
library(stringr)
library(yaml)

get_script <- function(filename) {
    lines <- readLines(filename, warn=FALSE)
    lines <- str_trim(lines, "right")

    state <- "main"
    need_blank <- FALSE
    output <- c()
    show <- function(...) output[length(output)+1] <<- paste0(...)
    for(line in lines) {
        if (state == "main" && str_detect(line,"^---")) {
            state <- "start"
        } else if (state == "start" && str_detect(line,"^---")) {
            state <- "main"
            need_blank <- FALSE
        } else if (state == "main" && str_detect(line,"^```.*include=FALSE")) {
            state <- "hidden_code"
        } else if (state == "main" && str_detect(line,"^```")) {
            state <- "code"
            show("")
        } else if (state %in% c("code","hidden_code") && str_detect(line,"^```")) {
            state <- "main"
            need_blank <- TRUE
        } else if (state == "main") {
            clean_line <- str_replace_all(line,"`|<details>|</details>|<summary>|</summary>|\\{\\..*\\}","")
            clean_line <- str_trim(clean_line, "right")
            if (clean_line=="" || clean_line=="\\" || clean_line=="***") {
                need_blank <- TRUE
            } else if (str_detect(clean_line,"^<!--.*-->$")) {
                # Comment.
            } else if (str_detect(clean_line,"^#+\\s+(\\w|\\()+")) { # Heading starts with a letter or round bracket
                # Heading
                clean_line <- str_replace_all(clean_line,"\\{.*\\}", "")
                clean_line <- str_trim(clean_line, "right")
                depth <- nchar(str_match(clean_line, "^#*")[1,1])
                if (depth == 1) clean_line <- paste0(clean_line, " ================")
                if (depth == 2) clean_line <- paste0(clean_line, " --------")
                show("")
                show("")
                show(clean_line)
                need_blank <- TRUE
            } else if (str_detect(clean_line,"^#+")) {
                # Blank heading
            } else {
                if (need_blank)
                    show("")
                show("# ", clean_line)
                need_blank <- FALSE
            }
        } else if (state == "code") {
            show(line)
        }
    }
    
    output
}


#filenames <- list.files(pattern="[0-9].*\\.Rmd") |> sort() |>
#    setdiff(c("01-01-setup.Rmd", "07-Acknowledgements.Rmd", "08-session_info.Rmd"))

filenames <- read_yaml("_bookdown.yml", readLines.warn=FALSE)$rmd_files |>
    setdiff(c(
        "index.Rmd", 
        "02-schedual.Rmd", 
        "01-01-setup.Rmd", 
        "05-01-resources.Rmd",
        "07-Acknowledgements.Rmd",
        "07-01-seuratobject.Rmd",
        "08-session_info.Rmd"))

output <- lapply(filenames, get_script) |> unlist()

writeLines(output, "docs/workshop.R")

# source("scripts/make_R_script.R")