# Learning State Verbs for TNA Simulation

Predefined collections of student learning action verbs organized by
category. These can be used as state names when generating TNA networks
to simulate realistic student learning behavior sequences.

## Usage

``` r
LEARNING_STATES
```

## Format

A named list with 8 categories of learning verbs.

## Details

The verbs are organized into seven categories based on learning science:

- metacognitive:

  Self-regulation and awareness actions: planning, monitoring,
  evaluating one's own learning.

- cognitive:

  Mental processing actions: understanding, analyzing, synthesizing
  information.

- behavioral:

  Observable study actions: reading, writing, practicing, note-taking.

- social:

  Interaction and collaboration actions: discussing, helping, teaching,
  asking.

- motivational:

  Effort and persistence actions: focusing, persisting, engaging,
  striving.

- affective:

  Emotional and attitudinal states: enjoying, coping, managing stress,
  curiosity.

- group_regulation:

  Socially shared regulation of learning (SSRL): collaborative planning,
  monitoring, adapting in groups.

- lms:

  Learning Management System actions: viewing content, accessing
  resources, submitting assignments, forum posts.

## Examples

``` r
# Get all available categories
names(get_learning_states())
#> NULL

# Get cognitive verbs only
get_learning_states("cognitive")
#>  [1] "Read"          "Study"         "Analyze"       "Summarize"    
#>  [5] "Memorize"      "Connect"       "Apply"         "Comprehend"   
#>  [9] "Synthesize"    "Compare"       "Contrast"      "Infer"        
#> [13] "Interpret"     "Elaborate"     "Encode"        "Retrieve"     
#> [17] "Process"       "Understand"    "Learn"         "Recognize"    
#> [21] "Recall"        "Integrate"     "Differentiate" "Abstract"     
#> [25] "Generalize"    "Classify"      "Categorize"    "Deduce"       
#> [29] "Reason"        "Conclude"     

# Get multiple categories
get_learning_states(c("metacognitive", "cognitive"))
#>  [1] "Plan"          "Monitor"       "Evaluate"      "Reflect"      
#>  [5] "Regulate"      "Adjust"        "Adapt"         "Check"        
#>  [9] "Assess"        "Judge"         "Strategize"    "Prioritize"   
#> [13] "Set_goals"     "Track"         "Self_assess"   "Calibrate"    
#> [17] "Diagnose"      "Forecast"      "Anticipate"    "Reconsider"   
#> [21] "Read"          "Study"         "Analyze"       "Summarize"    
#> [25] "Memorize"      "Connect"       "Apply"         "Comprehend"   
#> [29] "Synthesize"    "Compare"       "Contrast"      "Infer"        
#> [33] "Interpret"     "Elaborate"     "Encode"        "Retrieve"     
#> [37] "Process"       "Understand"    "Learn"         "Recognize"    
#> [41] "Recall"        "Integrate"     "Differentiate" "Abstract"     
#> [45] "Generalize"    "Classify"      "Categorize"    "Deduce"       
#> [49] "Reason"        "Conclude"     

# Get all verbs combined
get_learning_states("all")
#>   [1] "Plan"          "Monitor"       "Evaluate"      "Reflect"      
#>   [5] "Regulate"      "Adjust"        "Adapt"         "Check"        
#>   [9] "Assess"        "Judge"         "Strategize"    "Prioritize"   
#>  [13] "Set_goals"     "Track"         "Self_assess"   "Calibrate"    
#>  [17] "Diagnose"      "Forecast"      "Anticipate"    "Reconsider"   
#>  [21] "Read"          "Study"         "Analyze"       "Summarize"    
#>  [25] "Memorize"      "Connect"       "Apply"         "Comprehend"   
#>  [29] "Synthesize"    "Compare"       "Contrast"      "Infer"        
#>  [33] "Interpret"     "Elaborate"     "Encode"        "Retrieve"     
#>  [37] "Process"       "Understand"    "Learn"         "Recognize"    
#>  [41] "Recall"        "Integrate"     "Differentiate" "Abstract"     
#>  [45] "Generalize"    "Classify"      "Categorize"    "Deduce"       
#>  [49] "Reason"        "Conclude"      "Practice"      "Annotate"     
#>  [53] "Research"      "Review"        "Revise"        "Test"         
#>  [57] "Write"         "Note"          "Highlight"     "Underline"    
#>  [61] "Reread"        "Skim"          "Scan"          "Draft"        
#>  [65] "Edit"          "Copy"          "Record"        "Complete"     
#>  [69] "Submit"        "Attempt"       "Repeat"        "Drill"        
#>  [73] "Exercise"      "Rehearse"      "Outline"       "Diagram"      
#>  [77] "Map"           "List"          "Organize"      "Structure"    
#>  [81] "Collaborate"   "Discuss"       "Seek_help"     "Question"     
#>  [85] "Explain"       "Share"         "Teach"         "Tutor"        
#>  [89] "Debate"        "Argue"         "Negotiate"     "Consult"      
#>  [93] "Ask"           "Answer"        "Present"       "Participate"  
#>  [97] "Engage"        "Contribute"    "Support"       "Help"         
#> [101] "Clarify"       "Communicate"   "Listen"        "Respond"      
#> [105] "Feedback"      "Critique"      "Peer_review"   "Co_create"    
#> [109] "Brainstorm"    "Network"       "Focus"         "Persist"      
#> [113] "Explore"       "Create"        "Strive"        "Commit"       
#> [117] "Motivate"      "Endure"        "Overcome"      "Challenge"    
#> [121] "Aspire"        "Dedicate"      "Invest"        "Concentrate"  
#> [125] "Attend"        "Sustain"       "Maintain"      "Initiate"     
#> [129] "Continue"      "Pursue"        "Drive"         "Hustle"       
#> [133] "Push"          "Achieve"       "Accomplish"    "Excel"        
#> [137] "Improve"       "Grow"          "Develop"       "Progress"     
#> [141] "Enjoy"         "Appreciate"    "Value"         "Interest"     
#> [145] "Curious"       "Worry"         "Stress"        "Relax"        
#> [149] "Cope"          "Manage"        "Calm"          "Frustrate"    
#> [153] "Satisfy"       "Excite"        "Bore"          "Confuse"      
#> [157] "Resolve"       "Embrace"       "Accept"        "Tolerate"     
#> [161] "Celebrate"     "Doubt"         "Confident"     "Anxious"      
#> [165] "Hopeful"       "Discourage"    "Encourage"     "Inspire"      
#> [169] "Overwhelm"     "Relief"        "Cohesion"      "Consensus"    
#> [173] "Coregulate"    "Emotion"       "Synthesis"     "View"         
#> [177] "Access"        "Download"      "Upload"        "Click"        
#> [181] "Navigate"      "Browse"        "Login"         "Logout"       
#> [185] "Post"          "Reply"         "Forum"         "Quiz"         
#> [189] "Assignment"    "Video"         "Resource"      "Grade"        
#> [193] "Module"        "Page"          "File"          "Link"         
#> [197] "Course"        "Content"       "Discussion"    "Message"      
#> [201] "Announcement"  "Calendar"     

# Random selection of 8 verbs
get_learning_states(n = 8)
#> [1] "Calendar"      "Engage"        "Read"          "Differentiate"
#> [5] "Embrace"       "Reply"         "Prioritize"    "Submit"       

# Random 6 verbs from social and cognitive
get_learning_states(c("social", "cognitive"), n = 6)
#> [1] "Summarize"  "Generalize" "Infer"      "Co_create"  "Present"   
#> [6] "Seek_help" 
```
