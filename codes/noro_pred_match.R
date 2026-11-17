#  ──────────────────────────────────────────────────────────────────────────
#  Do the predicted mutations at a node turn up in the real tree?
#
#  Input
#    ~/noro_nodeall_long_final_unique.rds   129 nodes, one row each
#      node             integer node id
#      id               representative accession
#      date             decimal year
#      clinical_state   acute / chronic / unknown
#      id_list          descendant accessions (2-1182)
#      mutation_real_u  union of observed mutations over those descendants
#      mutation_pred    13-38 predicted mutations
#    offspring, pos, label and mutation_real are ignored.
#
#  Usage
#    source("noro_pred_match.R")
#    res <- pred_match()        # analysis, CSVs, and the numbers printed
#    figs <- pred_match_figs(res)   # four figures, saved as PNG + PDF
#
#  ── THE MATCHING RULE, AND WHY IT IS NOT STRING EQUALITY ──────────────────
#
#  A mutation is written <ancestral><position><derived>, e.g. T294P. The two
#  sets are anchored differently: observed calls are relative to the RefSeq,
#  predicted calls relative to the node's own reconstructed sequence. The
#  ancestral residue therefore disagrees whenever the node has drifted from
#  RefSeq at that codon, and T294P vs I294P would be scored a miss although
#  both say "position 294 becomes P". So the ancestral letter is dropped and
#  matching is done two ways:
#
#    position           the right codon was named
#    position+residue   the right codon AND the right destination residue
#
#  position+residue is a strict subset of position, which is what lets the
#  figures stack them as a partition of each node's predictions:
#
#    predictions = [position+residue] + [position only] + [no match]
#
#  Each node's predictions are de-duplicated to unique (position, residue)
#  pairs first -- otherwise a node predicting the same substitution twice
#  would score it twice.
#
#  ── THE TWO COMPARISON SETS ───────────────────────────────────────────────
#
#    own descendants   that node's own mutation_real_u
#    later nodes       union of mutation_real_u over every node with a
#                      STRICTLY greater date; the node itself is excluded
#
#  Both percentages are out of the node's own predictions -- precision, not
#  recall. "Of what this node predicted, how much showed up?"
#
#  A warning about the later-node set, which the by-date figure is there to
#  make visible: it shrinks monotonically as the node gets younger, and is
#  empty for the most recent node. A high match rate at an old node is partly
#  just a larger set to match against. Read the two panels of Figure 3
#  together, never the top one alone.
#
#  ── DOMAINS ───────────────────────────────────────────────────────────────
#
#  VP1 is 540 codons: S 1-225, P1 226-278 and 406-540, P2 279-405,
#  P domain = P1 + P2. The restricted analysis keeps predictions whose
#  position falls in the P domain, and the denominator is those predictions.
#  ──────────────────────────────────────────────────────────────────────────

library(tidyverse)


# ── CONFIG ────────────────────────────────────────────────────────────────

PM_IN  <- "~/noro_nodeall_long_final_unique.rds"
PM_DIR <- "/Users/rubayetalam/claudetry"
PM_STEM <- "noro_pred_match"

VP1_LEN <- 540
PM_DOMAINS <- tribble(
  ~domain, ~start, ~end,
  "S",       1,    225,
  "P1",    226,    278,
  "P2",    279,    405,
  "P1",    406,    540
)
P_DOMAIN <- c(226, 540)        # P1 + P2 as one contiguous span

#  Antigenic sites, same definition as noro_mutation.R.
PM_SITES <- list(
  A = c(294:298, 368, 372, 373), B = c(333, 389),
  C = c(339:341, 375:378),       D = c(393:398),
  E = c(407:414),                G = c(352, 355:357, 359, 364)
)

#  Three outcomes that partition every prediction, ordered worst to best so
#  the stack reads upward. One hue family, light to dark, because this is an
#  ordinal scale, not six unrelated categories.
PM_PAL <- c("No match"           = "#d9d8d4",
            "Position only"      = "#8fc7ea",
            "Position + residue" = "#1b5e9e")

#  Clinical state, identical to noro_phylo.R so the colours mean the same
#  thing across every figure in the project.
PM_TRAIT <- c(acute = "#2166ac", chronic = "#d6604d",
              persistent = "#1a9850", unknown = "grey30")

PM_INK <- c(primary = "#0b0b0b", secondary = "#3d3c3a",
            muted = "#7a7974", rule = "#2b2a28")

PM_DPI <- 600


# ── PARSING ───────────────────────────────────────────────────────────────

MUT_RX <- "^([A-Z])([0-9]+)([A-Z*])$"

#' Split mutation strings into position and derived residue.
#'
#' Anything that does not match MUT_RX is returned in `bad` rather than
#' silently dropped, so a malformed call is reported instead of quietly
#' lowering a denominator.
parse_mut <- function(v) {
  v <- v[!is.na(v) & nzchar(v)]
  ok <- grepl(MUT_RX, v)
  list(
    tbl = tibble(pos = as.integer(sub(MUT_RX, "\\2", v[ok])),
                 der = sub(MUT_RX, "\\3", v[ok]),
                 key = paste0(sub(MUT_RX, "\\2", v[ok]),
                              sub(MUT_RX, "\\3", v[ok]))),
    bad = v[!ok]
  )
}

in_p_domain <- function(pos) pos >= P_DOMAIN[1] & pos <= P_DOMAIN[2]

site_of <- function(pos, sites = PM_SITES) {
  out <- rep(NA_character_, length(pos))
  for (s in names(sites)) out[pos %in% sites[[s]]] <- s
  out
}


# ── ANALYSIS ──────────────────────────────────────────────────────────────

#' Count the two kinds of match of `p` (a node's predictions) in `s` (a
#' comparison set), and return the matched substitutions as a string.
#'
#' `s` may be empty -- the most recent node has no later nodes at all. That
#' is zero matches out of a set of size zero, not NA: the prediction was made
#' and nothing confirmed it.
match_one <- function(p, s) {
  if (!nrow(p))
    return(tibble(n_pred = 0L, pos_hit = 0L, both_hit = 0L,
                  matched = NA_character_))
  hit_pos  <- p$pos %in% s$pos
  hit_both <- p$key %in% s$key
  tibble(
    n_pred   = nrow(p),
    pos_hit  = sum(hit_pos),
    both_hit = sum(hit_both),
    #  The substitutions that matched on position AND residue, written back
    #  in the usual form minus the (incomparable) ancestral letter.
    matched  = if (any(hit_both))
                 paste(paste0(p$pos[hit_both], p$der[hit_both]), collapse = ";")
               else NA_character_
  )
}

#' Everything: parse, build both comparison sets, score, write the CSVs.
pred_match <- function(path = PM_IN, write = TRUE) {

  y <- readRDS(path.expand(path))

  pred_p <- map(y$mutation_pred,   parse_mut)
  real_p <- map(y$mutation_real_u, parse_mut)

  #  Flag parse failures loudly rather than letting them vanish.
  bad <- tibble(
    node = rep(y$node, lengths(map(pred_p, "bad")) + lengths(map(real_p, "bad"))),
    value = c(unlist(map(pred_p, "bad")), unlist(map(real_p, "bad")))
  )
  if (nrow(bad)) {
    warning(nrow(bad), " mutation string(s) did not parse; see $unparsed",
            call. = FALSE)
    print(bad)
  } else {
    message("all mutation strings parsed cleanly")
  }

  #  De-duplicate predictions to unique (position, residue) pairs.
  pred <- map(pred_p, ~ distinct(.x$tbl, pos, der, .keep_all = TRUE))
  real <- map(real_p, ~ distinct(.x$tbl, pos, der, .keep_all = TRUE))

  #  Comparison set 2: everything observed anywhere later in time. Strictly
  #  greater, so a node never matches against itself, and ties in date
  #  (same reconstructed year) are excluded as well -- they are not "later".
  later <- map(y$date, function(d) {
    i <- which(y$date > d)
    if (!length(i)) real[[1]][0, ] else distinct(bind_rows(real[i]), pos, der,
                                                 .keep_all = TRUE)
  })

  sets <- list(`Own descendants` = real, `Later nodes` = later)

  #  One row per node per comparison set, with the whole-VP1 figures and the
  #  P-domain-only figures side by side.
  rows <- imap_dfr(sets, function(set_list, set_name) {
    pmap_dfr(list(seq_len(nrow(y)), pred, set_list), function(i, p, s) {
      all_m <- match_one(p, s)
      pP    <- p[in_p_domain(p$pos), ]
      p_m   <- match_one(pP, s)
      tibble(
        node           = y$node[i],
        id             = y$id[i],
        date           = y$date[i],
        clinical_state = y$clinical_state[i],
        n_descendants  = length(y$id_list[[i]]),
        comparison     = set_name,
        n_set          = nrow(s),
        n_pred         = all_m$n_pred,
        pos_match      = all_m$pos_hit,
        pos_match_pct  = 100 * all_m$pos_hit / all_m$n_pred,
        posres_match   = all_m$both_hit,
        posres_match_pct = 100 * all_m$both_hit / all_m$n_pred,
        n_pred_P       = p_m$n_pred,
        pos_match_P    = p_m$pos_hit,
        pos_match_P_pct = 100 * p_m$pos_hit / p_m$n_pred,
        posres_match_P = p_m$both_hit,
        posres_match_P_pct = 100 * p_m$both_hit / p_m$n_pred,
        matched_subs   = all_m$matched
      )
    })
  }) %>%
    mutate(comparison = factor(comparison, levels = names(sets))) %>%
    arrange(comparison, date)

  #  ---- totals -----------------------------------------------------------
  #  Pooled over nodes, i.e. a prediction-weighted rate, not the mean of the
  #  per-node percentages. The per-node mean would give a node with 13
  #  predictions the same weight as one with 38.
  totals <- bind_rows(
    rows %>% group_by(comparison) %>%
      summarise(scope = "Whole VP1",
                n_pred = sum(n_pred), pos = sum(pos_match),
                posres = sum(posres_match), .groups = "drop"),
    rows %>% group_by(comparison) %>%
      summarise(scope = "P1 + P2",
                n_pred = sum(n_pred_P), pos = sum(pos_match_P),
                posres = sum(posres_match_P), .groups = "drop")
  ) %>%
    mutate(pos_pct = 100 * pos / n_pred, posres_pct = 100 * posres / n_pred,
           scope = factor(scope, levels = c("Whole VP1", "P1 + P2"))) %>%
    arrange(scope, comparison)

  by_state <- rows %>%
    group_by(comparison, clinical_state) %>%
    summarise(n_nodes = n(), n_pred = sum(n_pred),
              pos = sum(pos_match), posres = sum(posres_match),
              n_pred_P = sum(n_pred_P), pos_P = sum(pos_match_P),
              posres_P = sum(posres_match_P), .groups = "drop") %>%
    mutate(pos_pct = 100 * pos / n_pred, posres_pct = 100 * posres / n_pred,
           pos_P_pct = 100 * pos_P / n_pred_P,
           posres_P_pct = 100 * posres_P / n_pred_P)

  #  Per-position detail, for the protein map.
  by_pos <- imap_dfr(sets, function(set_list, set_name) {
    pmap_dfr(list(pred, set_list), function(p, s) {
      if (!nrow(p)) return(tibble())
      p %>% mutate(comparison = set_name,
                   outcome = case_when(key %in% s$key ~ "Position + residue",
                                       pos %in% s$pos ~ "Position only",
                                       TRUE           ~ "No match"))
    })
  }) %>%
    mutate(comparison = factor(comparison, levels = names(sets)),
           outcome = factor(outcome, levels = names(PM_PAL)),
           site = site_of(pos))

  out <- list(rows = rows, totals = totals, by_state = by_state,
              by_pos = by_pos, unparsed = bad)

  if (write) {
    write_csv(rows,     file.path(PM_DIR, paste0(PM_STEM, "_by_node.csv")))
    write_csv(totals,   file.path(PM_DIR, paste0(PM_STEM, "_totals.csv")))
    write_csv(by_state, file.path(PM_DIR, paste0(PM_STEM, "_by_state.csv")))
    write_csv(by_pos,   file.path(PM_DIR, paste0(PM_STEM, "_by_position.csv")))
    message("wrote ", PM_STEM, "_{by_node,totals,by_state,by_position}.csv")
  }

  print_pred_match(out)
  invisible(out)
}

print_pred_match <- function(res) {
  cat("\n== Totals (percentages are of predictions: precision) ==\n")
  print(res$totals %>%
          transmute(scope, comparison, n_pred,
                    position = sprintf("%d (%.1f%%)", pos, pos_pct),
                    `position + residue` = sprintf("%d (%.1f%%)", posres,
                                                   posres_pct)) %>%
          as.data.frame(), row.names = FALSE)

  cat("\n== By clinical state ==\n")
  print(res$by_state %>%
          transmute(comparison, clinical_state, n_nodes, n_pred,
                    position = sprintf("%d (%.1f%%)", pos, pos_pct),
                    `position + residue` = sprintf("%d (%.1f%%)", posres,
                                                   posres_pct)) %>%
          as.data.frame(), row.names = FALSE)
  invisible(res)
}


# ── THEME ─────────────────────────────────────────────────────────────────

theme_pm <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text          = element_text(colour = PM_INK[["primary"]]),
      plot.title    = element_text(face = "bold", size = rel(1.2),
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(colour = PM_INK[["secondary"]],
                                   size = rel(0.95), margin = margin(b = 8)),
      plot.caption  = element_text(colour = PM_INK[["muted"]], size = rel(0.8),
                                   hjust = 0, margin = margin(t = 8)),
      axis.line     = element_line(colour = PM_INK[["rule"]], linewidth = 0.3),
      axis.ticks    = element_line(colour = PM_INK[["rule"]], linewidth = 0.3),
      axis.text     = element_text(colour = PM_INK[["secondary"]]),
      panel.grid    = element_blank(),
      strip.background = element_blank(),
      strip.text    = element_text(face = "bold", size = rel(1)),
      legend.position = "top",
      legend.title  = element_blank(),
      legend.key.size = unit(10, "pt"),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}


# ── FIGURES ───────────────────────────────────────────────────────────────

#  Figure 1. The headline: what share of each node's predictions turns up in
#  each comparison set, at each matching strictness.
#
#  A stacked 100% bar was the obvious form and is the wrong one: three of the
#  four bars are over 90% "no match", so the bar spends its height on the
#  uninformative part and crushes the 1.3% and 1.6% slivers into nothing. A
#  paired dot on a percentage axis gives every value the same room. The two
#  dots on a row are joined because position+residue is nested inside
#  position -- the segment IS the cost of demanding the right residue too.
fig_overall <- function(res) {
  d <- res$totals %>%
    transmute(scope, comparison, n_pred,
              Position = pos, `Position + residue` = posres,
              p1 = pos_pct, p2 = posres_pct)

  seg <- d %>% transmute(scope, comparison, x = p2, xend = p1)

  long <- d %>%
    pivot_longer(c(Position, `Position + residue`),
                 names_to = "kind", values_to = "n") %>%
    mutate(pct = if_else(kind == "Position", p1, p2),
           kind = factor(kind, levels = c("Position", "Position + residue")))

  ggplot(long, aes(pct, comparison)) +
    geom_segment(data = seg, aes(x = x, xend = xend, y = comparison,
                                 yend = comparison),
                 inherit.aes = FALSE, linewidth = 0.9, colour = "#c7d9e8") +
    geom_point(aes(colour = kind), size = 3.6) +
    #  One label above the row, one below. Side by side they collide on the
    #  "Own descendants" rows, where the two values are four points apart.
    geom_text(aes(label = sprintf("%d  (%.1f%%)", n, pct), colour = kind,
                  vjust = if_else(kind == "Position", -1.3, 2.1)),
              size = 2.9, show.legend = FALSE) +
    facet_wrap(~ scope, ncol = 1) +
    scale_colour_manual(values = c(Position = "#8fc7ea",
                                   `Position + residue` = "#1b5e9e")) +
    scale_x_continuous(limits = c(-4, 58), expand = expansion(mult = c(0, 0.02)),
                       breaks = seq(0, 55, 10),
                       labels = scales::label_percent(scale = 1)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.75))) +
    labs(x = "Of this node\u2019s predictions", y = NULL,
         title = "Predicted mutations recovered in the real tree",
         subtitle = paste("2,770 unique predicted (position, residue) pairs",
                          "across 129 nodes; 1,178 of them in P1 + P2"),
         caption = paste0("Percentages are of predictions made (precision), ",
                          "not of mutations observed.\n",
                          "Position + residue is nested inside position: the ",
                          "connector is what demanding the right residue costs.")) +
    theme_pm() +
    theme(panel.grid.major.x = element_line(colour = "#eceae6",
                                            linewidth = 0.3),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

#  Figure 2. The same split, against later nodes only, by clinical state.
#  Wilson intervals on the position+residue rate, because n per state runs
#  from 101 to 2,166 predictions and the three rates are not equally certain.
fig_by_state <- function(res, conf = 0.95) {
  wilson <- function(k, n, conf = 0.95) {
    if (n == 0) return(c(NA, NA))
    z <- qnorm(1 - (1 - conf) / 2); p <- k / n
    c((p + z^2/(2*n) - z*sqrt((p*(1-p) + z^2/(4*n))/n)) / (1 + z^2/n),
      (p + z^2/(2*n) + z*sqrt((p*(1-p) + z^2/(4*n))/n)) / (1 + z^2/n))
  }

  d <- res$by_state %>% filter(comparison == "Later nodes") %>%
    mutate(lo = map2_dbl(posres, n_pred, ~ 100 * wilson(.x, .y, conf)[1]),
           hi = map2_dbl(posres, n_pred, ~ 100 * wilson(.x, .y, conf)[2]),
           lo_pos = map2_dbl(pos, n_pred, ~ 100 * wilson(.x, .y, conf)[1]),
           hi_pos = map2_dbl(pos, n_pred, ~ 100 * wilson(.x, .y, conf)[2]),
           clinical_state = fct_reorder(clinical_state, -posres_pct))

  long <- d %>%
    transmute(clinical_state, n_nodes, n_pred,
              Position = pos_pct, `Position + residue` = posres_pct,
              l1 = lo_pos, h1 = hi_pos, l2 = lo, h2 = hi) %>%
    pivot_longer(c(Position, `Position + residue`),
                 names_to = "kind", values_to = "pct") %>%
    mutate(lo = if_else(kind == "Position", l1, l2),
           hi = if_else(kind == "Position", h1, h2),
           kind = factor(kind, levels = c("Position", "Position + residue")))

  ggplot(long, aes(clinical_state, pct, fill = kind)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62,
             colour = "white", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.72), width = 0.14,
                  linewidth = 0.4, colour = PM_INK[["rule"]]) +
    geom_text(aes(y = hi, label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.72), vjust = -0.6,
              size = 3.0, colour = PM_INK[["secondary"]]) +
    geom_text(data = distinct(long, clinical_state, n_nodes, n_pred),
              aes(x = clinical_state, y = 0,
                  label = sprintf("%d nodes\n%d predictions", n_nodes, n_pred)),
              inherit.aes = FALSE, vjust = 1.4, size = 2.7,
              colour = PM_INK[["muted"]]) +
    scale_fill_manual(values = c(Position = "#8fc7ea",
                                 `Position + residue` = "#1b5e9e")) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.10)),
                       labels = scales::label_percent(scale = 1)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Predictions recovered later",
         title = "Recovery against later nodes, by clinical state",
         subtitle = sprintf("Whiskers: %.0f%% Wilson interval", 100 * conf),
         caption = paste0("'unknown' is 5 nodes and 101 predictions -- the ",
                          "interval is wide for a reason.\nState is the ",
                          "node's reconstructed state, not the state of the ",
                          "sequences it is matched against.")) +
    theme_pm()
}

#  Figure 3. Per-node rate against time, over the size of the comparison set.
#
#  Two panels sharing one x axis rather than one panel with two y axes: the
#  top rate is only interpretable next to the bottom count, and a dual axis
#  would invite reading a trend in the rate that is really a trend in how
#  much material was available to match.
fig_by_date <- function(res) {
  d <- res$rows %>% filter(comparison == "Later nodes")

  top <- ggplot(d, aes(date, posres_match_pct)) +
    geom_point(aes(size = n_pred, colour = clinical_state), alpha = 0.85) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = PM_INK[["rule"]], linewidth = 0.5, span = 0.8) +
    scale_colour_manual(values = PM_TRAIT) +
    scale_size_area(max_size = 4.5, breaks = c(15, 25, 35),
                    name = "predictions") +
    scale_y_continuous(labels = scales::label_percent(scale = 1),
                       expand = expansion(mult = c(0.03, 0.08))) +
    #  The loess runs below zero at the right-hand edge, where the later-node
    #  set empties out. A negative recovery rate is not a thing; clip it.
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = NULL, y = "Position + residue recovered",
         title = "A node's recovery rate tracks how much tree is left after it",
         subtitle = paste("Top: share of each node's predictions seen in any",
                          "later node. Bottom: how many distinct mutations",
                          "that set contains.")) +
    theme_pm() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank())

  bottom <- ggplot(d, aes(date, n_set)) +
    geom_area(fill = "#d9d8d4", alpha = 0.7) +
    geom_line(colour = PM_INK[["rule"]], linewidth = 0.4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Node date (decimal year)",
         y = "Mutations in\nthe later set",
         caption = paste0("The later-node set shrinks to nothing at the most ",
                          "recent node, so a low rate on the right is a ",
                          "smaller target, not a worse prediction.")) +
    theme_pm()

  patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(2.2, 1))
}

#  Figure 4. Where along VP1 the predictions land, and which of them are
#  recovered. Counts, not rates, so the eye is not drawn to a position with
#  one prediction and one match.
fig_position <- function(res) {
  d <- res$by_pos %>% filter(comparison == "Later nodes") %>%
    count(pos, outcome)

  dom <- PM_DOMAINS %>% mutate(mid = (start + end) / 2)
  sites <- imap_dfr(PM_SITES, ~ tibble(site = .y, pos = .x))
  ymax <- d %>% group_by(pos) %>% summarise(n = sum(n)) %>% pull(n) %>% max()

  ggplot() +
    #  Domain bands behind everything, labelled once each.
    geom_rect(data = dom, aes(xmin = start - 0.5, xmax = end + 0.5,
                              ymin = -Inf, ymax = Inf, alpha = domain),
              fill = "#9aa0a6") +
    scale_alpha_manual(values = c(S = 0.05, P1 = 0.11, P2 = 0.20),
                       guide = "none") +
    #  Antigenic-site residues as ticks along the baseline.
    geom_point(data = sites, aes(x = pos, y = -ymax * 0.045),
               shape = 15, size = 0.9, colour = "#eb6834") +
    geom_col(data = d, aes(pos, n, fill = outcome), width = 1.1) +
    scale_fill_manual(values = PM_PAL, breaks = rev(names(PM_PAL))) +
    annotate("text", x = dom$mid, y = ymax * 1.06, label = dom$domain,
             fontface = "bold", size = 3.1, colour = PM_INK[["secondary"]]) +
    scale_x_continuous(breaks = seq(0, VP1_LEN, 60), limits = c(0, VP1_LEN + 1),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "VP1 codon position", y = "Predictions at this position",
         title = "Where the predictions are made, and where they are recovered",
         subtitle = paste("All 2,770 predictions, scored against later nodes.",
                          "Orange ticks mark the 37 blockade-epitope residues."),
         caption = paste0("Predictions are densest in S -- 21.2 per residue ",
                          "against 9.2 in P2 -- but recover worst there ",
                          "(14.9% vs 25.7%).\nOn the 37 epitope residues ",
                          "52.3% of predictions recover at position+residue, ",
                          "against 14.2% elsewhere.")) +
    theme_pm() +
    guides(fill = guide_legend(reverse = TRUE))
}

#' Build all four figures and save each as PNG (600 dpi) and PDF.
pred_match_figs <- function(res, save = TRUE) {
  figs <- list(
    fig1_overall  = list(p = fig_overall(res),  w = 7.6, h = 5.2),
    fig2_by_state = list(p = fig_by_state(res), w = 7.0, h = 5.4),
    fig3_by_date  = list(p = fig_by_date(res),  w = 8.2, h = 6.4),
    fig4_position = list(p = fig_position(res), w = 10.5, h = 4.8)
  )

  if (save) {
    dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
    for (nm in names(figs)) {
      f <- figs[[nm]]
      stem <- file.path(PM_DIR, paste0(PM_STEM, "_", nm))
      ggsave(paste0(stem, ".png"), f$p, device = dev, width = f$w, height = f$h,
             dpi = PM_DPI, bg = "white")
      ggsave(paste0(stem, ".pdf"), f$p, device = cairo_pdf, width = f$w,
             height = f$h, bg = "white")
    }
    message("wrote ", length(figs), " figures as .png (", PM_DPI,
            " dpi) and .pdf")
  }
  invisible(map(figs, "p"))
}


# ── EXAMPLES ──────────────────────────────────────────────────────────────
#
# source("noro_pred_match.R")
# res  <- pred_match()            # CSVs + the printed tables
# figs <- pred_match_figs(res)    # four figures
# figs$fig1_overall               # redraw one on screen
#
# The per-node table is the thing to join onto other analyses:
# res$rows %>% filter(comparison == "Later nodes") %>%
#   arrange(desc(posres_match_pct)) %>%
#   select(node, date, clinical_state, n_pred, posres_match, matched_subs)
