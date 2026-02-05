#' Learning State Verbs for TNA Simulation
#'
#' @description
#' Predefined collections of student learning action verbs organized by
#' category. These can be used as state names when generating TNA networks
#' to simulate realistic student learning behavior sequences.
#'
#' @details
#' The verbs are organized into seven categories based on learning science:
#'
#' \describe{
#'   \item{metacognitive}{Self-regulation and awareness actions:
#'     planning, monitoring, evaluating one's own learning.}
#'   \item{cognitive}{Mental processing actions:
#'     understanding, analyzing, synthesizing information.}
#'   \item{behavioral}{Observable study actions:
#'     reading, writing, practicing, note-taking.}
#'   \item{social}{Interaction and collaboration actions:
#'     discussing, helping, teaching, asking.}
#'   \item{motivational}{Effort and persistence actions:
#'     focusing, persisting, engaging, striving.}
#'   \item{affective}{Emotional and attitudinal states:
#'     enjoying, coping, managing stress, curiosity.}
#'   \item{group_regulation}{Socially shared regulation of learning (SSRL):
#'     collaborative planning, monitoring, adapting in groups.}
#'   \item{lms}{Learning Management System actions:
#'     viewing content, accessing resources, submitting assignments, forum posts.}
#' }
#'
#' @examples
#' # Get all available categories
#' names(get_learning_states())
#'
#' # Get cognitive verbs only
#' get_learning_states("cognitive")
#'
#' # Get multiple categories
#' get_learning_states(c("metacognitive", "cognitive"))
#'
#' # Get all verbs combined
#' get_learning_states("all")
#'
#' # Random selection of 8 verbs
#' get_learning_states(n = 8)
#'
#' # Random 6 verbs from social and cognitive
#' get_learning_states(c("social", "cognitive"), n = 6)
#'
#' @name learning_states
#' @rdname learning_states
NULL

#' @rdname learning_states
#' @format A named list with 8 categories of learning verbs.
#' @export
LEARNING_STATES <- list(
  metacognitive = c(
    "Plan",
    "Monitor",
    "Evaluate",
    "Reflect",
    "Regulate",
    "Adjust",
    "Adapt",
    "Check",
    "Assess",
    "Judge",
    "Strategize",
    "Prioritize",
    "Set_goals",
    "Track",
    "Self_assess",
    "Calibrate",
    "Diagnose",
    "Forecast",
    "Anticipate",
    "Reconsider"
  ),

  cognitive = c(
    "Read",
    "Study",
    "Analyze",
    "Summarize",
    "Memorize",
    "Connect",
    "Apply",
    "Comprehend",
    "Synthesize",
    "Compare",
    "Contrast",
    "Infer",
    "Interpret",
    "Elaborate",
    "Encode",
    "Retrieve",
    "Process",
    "Understand",
    "Learn",
    "Recognize",
    "Recall",
    "Integrate",
    "Differentiate",
    "Abstract",
    "Generalize",
    "Classify",
    "Categorize",
    "Deduce",
    "Reason",
    "Conclude"
  ),

  behavioral = c(
    "Practice",
    "Annotate",
    "Research",
    "Review",
    "Revise",
    "Test",
    "Write",
    "Note",
    "Highlight",
    "Underline",
    "Reread",
    "Skim",
    "Scan",
    "Draft",
    "Edit",
    "Copy",
    "Record",
    "Complete",
    "Submit",
    "Attempt",
    "Repeat",
    "Drill",
    "Exercise",
    "Rehearse",
    "Outline",
    "Diagram",
    "Map",
    "List",
    "Organize",
    "Structure"
  ),

  social = c(
    "Collaborate",
    "Discuss",
    "Seek_help",
    "Question",
    "Explain",
    "Share",
    "Teach",
    "Tutor",
    "Debate",
    "Argue",
    "Negotiate",
    "Consult",
    "Ask",
    "Answer",
    "Present",
    "Participate",
    "Engage",
    "Contribute",
    "Support",
    "Help",
    "Clarify",
    "Communicate",
    "Listen",
    "Respond",
    "Feedback",
    "Critique",
    "Peer_review",
    "Co_create",
    "Brainstorm",
    "Network"
  ),

  motivational = c(
    "Focus",
    "Persist",
    "Explore",
    "Create",
    "Strive",
    "Commit",
    "Motivate",
    "Endure",
    "Overcome",
    "Challenge",
    "Aspire",
    "Dedicate",
    "Invest",
    "Concentrate",
    "Attend",
    "Sustain",
    "Maintain",
    "Initiate",
    "Continue",
    "Pursue",
    "Drive",
    "Hustle",
    "Push",
    "Achieve",
    "Accomplish",
    "Excel",
    "Improve",
    "Grow",
    "Develop",
    "Progress"
  ),

  affective = c(
    "Enjoy",
    "Appreciate",
    "Value",
    "Interest",
    "Curious",
    "Worry",
    "Stress",
    "Relax",
    "Cope",
    "Manage",
    "Calm",
    "Frustrate",
    "Satisfy",
    "Excite",
    "Bore",
    "Confuse",
    "Resolve",
    "Embrace",
    "Accept",
    "Tolerate",
    "Celebrate",
    "Doubt",
    "Confident",
    "Anxious",
    "Hopeful",
    "Discourage",
    "Encourage",
    "Inspire",
    "Overwhelm",
    "Relief"
  ),

  group_regulation = c(
    "Adapt",
    "Cohesion",
    "Consensus",
    "Coregulate",
    "Discuss",
    "Emotion",
    "Monitor",
    "Plan",
    "Synthesis"
  ),

  lms = c(
    "View",
    "Access",
    "Download",
    "Upload",
    "Submit",
    "Click",
    "Navigate",
    "Browse",
    "Login",
    "Logout",
    "Post",
    "Reply",
    "Forum",
    "Quiz",
    "Assignment",
    "Video",
    "Resource",
    "Grade",
    "Attempt",
    "Complete",
    "Module",
    "Page",
    "File",
    "Link",
    "Course",
    "Content",
    "Discussion",
    "Message",
    "Announcement",
    "Calendar"
  )
)


#' Get Learning State Verbs
#'
#' @description
#' Retrieve learning action verbs by category, with options for random
#' selection and combining multiple categories.
#'
#' @param categories Character vector specifying which categories to include.
#'   Options: "metacognitive", "cognitive", "behavioral", "social",
#'   "motivational", "affective", or "all" for all categories.
#'   Default: "all".
#' @param n Integer or NULL. If specified, randomly sample n verbs from
#'   the selected categories. If NULL, return all verbs. Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducible sampling.
#'   Default: NULL.
#'
#' @return Character vector of learning state verbs.
#'
#' @examples
#' # All metacognitive verbs
#' get_learning_states("metacognitive")
#'
#' # 10 random verbs from all categories
#' get_learning_states(n = 10, seed = 42)
#'
#' # 5 verbs from cognitive and behavioral
#' get_learning_states(c("cognitive", "behavioral"), n = 5)
#'
#' # Self-regulated learning focus (metacognitive + motivational)
#' get_learning_states(c("metacognitive", "motivational"))
#'
#' @export
get_learning_states <- function(categories = "all", n = NULL, seed = NULL) {
  valid_categories <- c("metacognitive", "cognitive", "behavioral",
                        "social", "motivational", "affective",
                        "group_regulation", "lms", "all")

  # Handle "all" category

if ("all" %in% categories) {
    categories <- setdiff(names(LEARNING_STATES), "all")
  }

  # Validate categories
  invalid <- setdiff(categories, names(LEARNING_STATES))
  if (length(invalid) > 0) {
    stop(sprintf(
      "Invalid categories: %s. Valid options: %s",
      paste(invalid, collapse = ", "),
      paste(valid_categories, collapse = ", ")
    ))
  }

  # Combine verbs from selected categories
  verbs <- unlist(LEARNING_STATES[categories], use.names = FALSE)
  verbs <- unique(verbs)

  # Random sampling if n is specified
  if (!is.null(n)) {
    if (!is.null(seed)) set.seed(seed)
    if (n > length(verbs)) {
      warning(sprintf(
        "Requested %d verbs but only %d available. Returning all.",
        n, length(verbs)
      ))
      n <- length(verbs)
    }
    verbs <- sample(verbs, n)
  }

  return(verbs)
}


#' Get Learning State Categories Summary
#'
#' @description
#' Display a summary of all available learning state categories and
#' the number of verbs in each.
#'
#' @return A data frame with category names and verb counts.
#'
#' @examples
#' list_learning_categories()
#'
#' @export
list_learning_categories <- function() {
  categories <- names(LEARNING_STATES)
  counts <- sapply(LEARNING_STATES, length)
  examples <- sapply(LEARNING_STATES, function(x) {
    paste(head(x, 5), collapse = ", ")
  })

  df <- data.frame(
    Category = categories,
    Count = counts,
    Examples = paste0(examples, ", ..."),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(df)
}


#' Smart State Selection for Networks
#'
#' @description
#' Intelligently select learning states based on the number of nodes needed,
#' optionally biasing toward specific categories.
#'
#' @param n_states Integer. Number of states needed.
#' @param primary_categories Character vector. Categories to prioritize.
#'   Default: NULL (balanced selection).
#' @param secondary_categories Character vector. Categories to supplement with.
#'   Default: NULL (use all remaining).
#' @param primary_ratio Numeric (0 to 1). Proportion of states from primary
#'   categories. Default: 0.6.
#' @param seed Integer or NULL. Random seed. Default: NULL.
#'
#' @return Character vector of selected learning states.
#'
#' @details
#' Selection logic:
#' \itemize{
#'   \item If n_states <= 5: Single category or balanced small set
#'   \item If n_states 6-10: 1-2 categories prioritized
#'   \item If n_states 11-20: Multiple categories with primary focus
#'   \item If n_states > 20: All categories combined
#' }
#'
#' @examples
#' # 5 states focused on self-regulation
#' smart_select_states(5, primary_categories = "metacognitive")
#'
#' # 10 states: mostly cognitive, some behavioral
#' smart_select_states(10,
#'   primary_categories = "cognitive",
#'   secondary_categories = "behavioral"
#' )
#'
#' # 15 states: balanced across categories
#' smart_select_states(15, seed = 42)
#'
#' @export
smart_select_states <- function(n_states,
                                 primary_categories = NULL,
                                 secondary_categories = NULL,
                                 primary_ratio = 0.6,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  all_categories <- names(LEARNING_STATES)

  # Auto-select strategy based on n_states if no categories specified
  if (is.null(primary_categories)) {
    if (n_states <= 5) {
      # Small network: pick one category
      primary_categories <- sample(all_categories, 1)
    } else if (n_states <= 10) {
      # Medium network: pick 2 categories
      primary_categories <- sample(all_categories, 2)
    } else if (n_states <= 20) {
      # Large network: pick 3 categories
      primary_categories <- sample(all_categories, 3)
    } else {
      # Very large: use all
      primary_categories <- all_categories
    }
  }

  # Determine secondary categories
  if (is.null(secondary_categories)) {
    secondary_categories <- setdiff(all_categories, primary_categories)
  }

  # Calculate how many from each group
  n_primary <- ceiling(n_states * primary_ratio)
  n_secondary <- n_states - n_primary

  # Get verbs from each group
  primary_verbs <- get_learning_states(primary_categories)
  secondary_verbs <- get_learning_states(secondary_categories)

  # Adjust if not enough verbs
  n_primary <- min(n_primary, length(primary_verbs))
  n_secondary <- min(n_secondary, length(secondary_verbs))

  # If still not enough, take more from primary
  if (n_primary + n_secondary < n_states) {
    n_primary <- min(n_states - n_secondary, length(primary_verbs))
  }

  # Sample
  selected_primary <- sample(primary_verbs, n_primary)
  selected_secondary <- if (n_secondary > 0 && length(secondary_verbs) > 0) {
    sample(secondary_verbs, min(n_secondary, length(secondary_verbs)))
  } else {
    character(0)
  }

  # Combine and shuffle
  selected <- c(selected_primary, selected_secondary)
  selected <- sample(selected, length(selected))

  # Ensure we have enough (pad with any remaining if needed)
  if (length(selected) < n_states) {
    all_verbs <- get_learning_states("all")
    remaining <- setdiff(all_verbs, selected)
    n_needed <- n_states - length(selected)
    if (length(remaining) >= n_needed) {
      selected <- c(selected, sample(remaining, n_needed))
    } else {
      selected <- c(selected, remaining)
      warning(sprintf(
        "Only %d unique verbs available, returning %d states.",
        length(selected), length(selected)
      ))
    }
  }

  return(selected[1:min(n_states, length(selected))])
}
