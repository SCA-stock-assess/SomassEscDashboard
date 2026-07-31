# ---- Shared helpers -------------------------------------------------------

ANCHOR_YEAR <- 2001  # fixed non-leap anchor used across all timing plots

# Leap-year-safe day-of-year: anchors every date to a fixed non-leap year
# before extracting day-of-year, so Aug-Oct dates align consistently across
# leap and non-leap calendar years. Replaces format(date, "%j").
safe_julian <- function(date) {
  lubridate::yday(as.Date(paste0(ANCHOR_YEAR, "-", format(date, "%m-%d"))))
}

# Converts a leap-safe julian day back to a plottable Date, anchored to the
# same fixed year used by safe_julian(). Shared by every timing plot so the
# x-axis is always built the same way.
julian_to_date <- function(j) as.Date(j - 1, origin = paste0(ANCHOR_YEAR, "-01-01"))

# Bounded [0,1] smoother used for the l95/u95 ribbon and the historic mean
# label line. Not a real binomial model (l95/u95 are proportions, not
# Bernoulli outcomes) -- used purely as a logit-constrained smoothing curve.
# The "non-integer #successes" warning this throws is expected and harmless
# for that reason; it's a signal you're using glm() outside its intended
# purpose, not a computation error.
logit_smooth <- function(y, x) {
  binomial()$linkinv(predict(glm(y ~ x, family = binomial)))
}

# Shared palette across all Stamp Chinook timing plots
hist_colour_dark  <- "#2C5F7C"  # muted steel blue (historic mean line / ribbon)
hist_colour_light <- "#A9C6D6"  # pale steel blue (unused now spaghetti uses a categorical palette)
curr_colour       <- "#2E8B57"  # sea green (current year -- reserved, never used for historic years)
target_colour     <- "#4D4D4D"  # dark grey (target line)

# Distinct, non-green palette for historic spaghetti-plot years. Green is
# deliberately excluded here so curr_colour above always reads unambiguously
# as "this year" against any number of historic lines.
hist_year_palette <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#B07AA1", "#9C755F",
  "#FF9DA7", "#BAB0AC", "#D4A6C8", "#FFBE7D", "#6B4C9A"
)

# Shared theme so every Stamp Chinook plot looks consistent
cn_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      axis.line.y = element_line(colour = "black", linewidth = 0.4),
      axis.line.x = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length.x = unit(0.15, "cm")
    )
}

# ---- Data loading -----------------------------------------------------

SomassEsc <- read_xlsx(
  "//dcbcpbsna01a.ENT.dfo-mpo.ca/PBS_SA_DFS$/SCD_Stad/WCVI/Sockeye/SOMASS/Data/ESCAPEMENT_PROGRAM/SomassEsc.xlsx",
  sheet = "Stamp CN & CO",
  na = ""
) |>
  select(1:9) |>
  pivot_longer(cols = Coho:ChinookJack, names_to = "species", values_to = "count") |>
  rename_with(tolower) |>
  mutate(
    date = as.Date(date),
    # Historic (complete) years: missing daily counts mean "no fish that
    # day", i.e. 0 -- not unknown. Leaves the current year's NAs alone
    # (handled separately via CurrentYearEsc / the live file).
    count = if_else(year < max(year) & is.na(count), 0, count)
  ) |>
  group_by(year, species) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    cum_count = cumsum(count),
    ann_ttl = sum(count, na.rm = TRUE),
    cum_prop = cum_count / ann_ttl,
    julian = safe_julian(date)
  ) |>
  ungroup()

CurrentYearEsc <- read_xlsx(
  "//dcbcpbsna01a.ENT.dfo-mpo.ca/PBS_SA_DFS$/SCD_Stad/WCVI/SOCKEYE/SOMASS/SOCKEYE_MGMT/2026_MGT/Daily Totals by Age 2026.xlsx",
  sheet = "Stamp CN&CO",
  na = ""
) |>
  select(1:9) |>
  pivot_longer(cols = 2:9, names_to = "species", values_to = "count") |>
  rename_with(tolower) |>
  filter(!is.na(date)) |>  # drop subtotal rows at end of sheet, if present
  mutate(date = as.Date(date)) |>
  group_by(species) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    cum_count = cumsum(count),
    ann_ttl = sum(count, na.rm = TRUE),
    cum_prop = cum_count / ann_ttl,
    julian = safe_julian(date)
  ) |>
  ungroup()

# ---- Timing (proportion) plot ------------------------------------------

#' Build Stamp River Chinook escapement timing plot
#'
#' @param hist_data    Historic data frame (e.g. SomassEsc). Columns needed:
#'                     species, year, julian, cum_prop.
#' @param current_data Current-season data frame (e.g. CurrentYearEsc),
#'                      updated weekly from the live file. Columns needed:
#'                      species, julian, cum_count.
#' @param curr_year    Current management year (integer), used for labels.
#' @param hist_years   Number of prior years in the historic mean/ribbon
#'                      (default 20).
#' @param esc_target   Numeric escapement target for curr_year.
build_cn_timing_plot <- function(hist_data, current_data, curr_year,
                                 hist_years = 20, esc_target) {
  
  hist_input <- hist_data |>
    filter(species == "Chinook",
           between(year, curr_year - hist_years, curr_year - 1))
  
  zero_total_years <- hist_input |>
    filter(!is.finite(cum_prop)) |>
    distinct(year) |>
    pull(year)
  
  if (length(zero_total_years) > 0) {
    warning(
      "Excluding year(s) with zero total Chinook count from historic average: ",
      paste(zero_total_years, collapse = ", ")
    )
  }
  
  hist_summary <- hist_input |>
    filter(is.finite(cum_prop)) |>
    group_by(julian) |>
    summarise(
      mean = mean(cum_prop),
      l95  = quantile(cum_prop, 0.05),
      u95  = quantile(cum_prop, 0.95),
      .groups = "drop"
    ) |>
    mutate(
      date = julian_to_date(julian),
      l95_smooth = logit_smooth(l95, julian),
      u95_smooth = logit_smooth(u95, julian)
    )
  
  current_input <- current_data |>
    filter(species == "Chinook") |>
    mutate(date = julian_to_date(julian), prop_of_target = cum_count / esc_target)
  
  if (nrow(current_input) == 0) {
    warning(
      "current_data has 0 rows after filtering species == 'Chinook' -- ",
      "check the actual species values with unique(current_data$species)."
    )
  }
  
  ggplot() +
    geom_ribbon(
      data = hist_summary,
      aes(date, ymin = l95_smooth, ymax = u95_smooth),
      fill = hist_colour_dark, alpha = 0.15
    ) +
    geom_textsmooth(
      data = hist_summary,
      aes(date, mean),
      label = paste("Historic", curr_year - hist_years, "\u2013", curr_year - 1),
      linewidth = 1, hjust = 0.6, colour = hist_colour_dark,
      method = "glm", method.args = list(family = binomial())
    ) +
    geom_hline(
      yintercept = 1, linetype = "dashed", colour = target_colour, linewidth = 0.5
    ) +
    geom_textline(
      data = current_input,
      aes(date, prop_of_target),
      label = as.character(curr_year),
      colour = curr_colour, hjust = 0.7, vjust = 0.8,
      linewidth = 1.1, text_smoothing = 60
    ) +
    scale_y_continuous(
      labels = scales::percent,
      name = "Proportion of total escapement (historic) / target (current year)",
      sec.axis = sec_axis(
        transform = ~ . * esc_target,
        labels = scales::comma,
        name = paste(curr_year, "cumulative escapement")
      ),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_date(breaks = "2 weeks", date_labels = "%d %b") +
    coord_cartesian(
      xlim = as.Date(c(paste0(ANCHOR_YEAR, "-08-01"), paste0(ANCHOR_YEAR, "-10-27")))
    ) +
    labs(x = NULL) +
    cn_theme() +
    theme(legend.position = "none")
}

# ---- Spaghetti (cumulative count) plot ----------------------------------

#' Build Stamp River Chinook cumulative-count spaghetti plot
#'
#' Shows each historic year as an individual labelled line (faded by
#' recency) plus the current year highlighted, on raw cumulative count
#' rather than proportion -- complements build_cn_timing_plot().
#'
#' @param hist_data    Historic data frame (e.g. SomassEsc). Columns needed:
#'                     species, year, julian, cum_count.
#' @param current_data Current-season data frame (e.g. CurrentYearEsc).
#'                      Columns needed: species, julian, cum_count.
#' @param curr_year    Current management year (integer), used for label.
#' @param hist_years   Number of prior years shown as individual lines
#'                      (default 10).
build_cn_spaghetti_plot <- function(hist_data, current_data, curr_year,
                                    hist_years = 10) {
  
  hist_input <- hist_data |>
    filter(species == "Chinook",
           between(year, curr_year - hist_years, curr_year - 1),
           julian < 310) |>
    mutate(date = julian_to_date(julian))
  
  current_input <- current_data |>
    filter(species == "Chinook", julian < 310) |>
    mutate(date = julian_to_date(julian))
  
  if (nrow(current_input) == 0) {
    warning(
      "current_data has 0 rows after filtering species == 'Chinook' -- ",
      "check the actual species values with unique(current_data$species)."
    )
  }
  
  hist_years_present <- sort(unique(hist_input$year))
  n_years <- length(hist_years_present)
  
  # Recycle/interpolate the palette to however many years are actually
  # present, so this works whether hist_years is 5 or 20.
  year_colours <- colorRampPalette(hist_year_palette)(n_years)
  names(year_colours) <- hist_years_present
  
  # Stagger each year's label along the line instead of clustering them all
  # near hjust = 0.95 -- that clustering was causing geomtextpath to drop
  # labels for lines that collided, which is why most years weren't showing.
  hist_input <- hist_input |>
    mutate(
      year_rank = match(year, hist_years_present),
      label_hjust = if (n_years > 1) {
        0.15 + 0.75 * (year_rank - 1) / (n_years - 1)
      } else {
        0.9
      }
    )
  
  ggplot(hist_input, aes(date, cum_count)) +
    geom_textline(
      aes(label = year, group = year, colour = factor(year), hjust = label_hjust),
      alpha = 0.85, linewidth = 0.6, size = 3, gap = FALSE
    ) +
    scale_colour_manual(values = year_colours, guide = "none") +
    geom_labelline(
      data = current_input,
      aes(date, cum_count),
      label = as.character(curr_year),
      colour = curr_colour,
      hjust = 0.65,
      linewidth = 1.5,
      boxcolour = "white",
      alpha = 0.9,
      label.padding = unit(0.1, "lines"),
      gap = TRUE
    ) +
    scale_x_date(
      breaks = "2 weeks", date_labels = "%d %b",
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(
      xlim = as.Date(c(paste0(ANCHOR_YEAR, "-08-01"), paste0(ANCHOR_YEAR, "-10-20")))
    ) +
    scale_y_continuous(
      position = "right",  # count values visible at the end of the time series
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = NULL, y = "Cumulative Stamp Falls Chinook escapement") +
    cn_theme() +
    theme(
      axis.title.y.right = element_text(margin = margin(l = 0.5, unit = "lines"))
    )
}

# ---- Running the functions and generating plots -----------
ChinookTimingPlot<- build_cn_timing_plot(
  hist_data    = SomassEsc,
  current_data = CurrentYearEsc,
  curr_year    = 2026,
  hist_years   = 20,
  esc_target   = 34000
)

ChinookSpagettiPlot<- build_cn_spaghetti_plot(
  hist_data    = SomassEsc,
  current_data = CurrentYearEsc,
  curr_year    = 2026,
  hist_years   = 10
)

ggsave(
  plot = ChinookSpagettiPlot,
  filename = paste0(
    "//dcbcpbsna01a.ENT.dfo-mpo.ca/PBS_SA_DFS$/SCD_Stad/WCVI/CHINOOK/CHINOOK_MGT/",
    curr_year,
    "/A23/Escapement plot/",
    "Fig4_2025_CN_Spagethi.png"
  ),
  height = 4.5,
  width = 8,
  units = "in"
)
