library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)

create_fixations_plot <- function(perm_results_signif, predictor_sep = ":", outfile = NULL) {
    fixations_plot <- perm_results_signif %>% 
        # make sure predictor is a character vector
        mutate(predictor = as.character(predictor)) %>%
        
        # compute the new level ordering:
        { 
            lvl <- unique(.$predictor)
            # count predictor sep characters in each label
            n_predictor_sep <- str_count(lvl, fixed(predictor_sep))
            # order by n_predictor_sep ascending, then by label alphabetically
            new_levels <- lvl[order(n_predictor_sep, lvl)]
            
            # then set predictor as a factor with those levels
            mutate(., predictor = factor(predictor, levels = new_levels))
        } %>%
        mutate(masked_p = ifelse(fdr_q_value > 0.05, NA, fdr_q_value)) %>%
        mutate(keep = case_when(gsub('baseline','',predictor) != predictor ~ FALSE,
                                TRUE ~ TRUE)) %>% # TODO comment out these two lines if you want to plot the baseline coefficient, too
        filter(keep == TRUE) %>% # TODO comment this line if you want to plot the baseline coefficient, too
        filter(time_bin >= 300 & time_bin <= 700) %>% #time_bins of interest: N4/P3: 300-500, later window 500-800
        drop_na() %>% 
        ggplot(., aes(x = time_bin, y = predictor, fill =fdr_q_value)) +
        geom_tile() +
        scale_fill_viridis(option = "viridis", name = "Permutation p-value (FDR corrected)")

    if (!is.null(outfile)) {
        ggsave(
            filename = outfile,
            plot = fixations_plot,
            width = 8,
            height = 6,
            dpi = 300
        )
    }

    invisible(fixations_plot)
}

