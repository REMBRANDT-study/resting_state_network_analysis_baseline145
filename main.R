# Analysis script for REMBRANDT network analysis with Higher Criticism
# Andrew Gerlach, 9/20/23

library(tidyverse)
library(readxl)
library(R.matlab)
library(DescTools)
library(reshape2)
library(oro.nifti)
library(survival)

# please set R working directory to the repository:
# setwd("<path>/resting_state_network_analysis_baseline145/")
source("calc_HC.R")
source("calc_HC_p_MC.R")

set.seed(123)

################################################################################
# Read in basic studay data (demographics, etc.)                               #
################################################################################

data_file <- "network_analysis_data.rds"
atlases <- c("Schaefer200", "Schaefer400", "Shen268")
networks <- c("Default", "Cont", "SalVentAttn", "DorsAttn", "Limbic", "Vis", "SomMot")
network_def_file <- "network_definitions/ATLAS_ROI_Network.txt"
test_types <- c("t_greater",
    "t_less",
    "reg_greater",
    "reg_less",
    "anova",
    "ancova",
    "reg_greater_nev_rem",
    "reg_less_nev_rem",
    "reg_greater_nev_rel",
    "reg_less_nev_rel",
    "reg_greater_rem_rel",
    "reg_less_rem_rel",
    "surv_greater",
    "surv_less")
n_atlas <- length(atlases)
n_net <- length(networks)
n_test <- length(test_types)

tmp <- readRDS(data_file)
dat <- tmp$dat
conn_matrices <- tmp$conn_matrices
n_subj <- nrow(dat)

################################################################################
# Analysis                                                                     #
################################################################################

p_table <- tibble(Atlas = character(),
    Network1 = character(),
    Network2 = character(),
    Node1 = numeric(),
    Node2 = numeric(),
    Test_type = character(),
    p = numeric())
hc_table <- tibble(Atlas = character(),
    Level = character(),
    Network1 = character(),
    Network2 = character(),
    Node = numeric(),
    Test_type = character(),
    hc = numeric(),
    p = numeric(),
    Bonferroni = numeric())

# perform analysis
for(atlas in 1:n_atlas) {

    cat(paste("Atlas:", atlases[atlas], "\n"))
    network_def <- read.csv(str_replace(network_def_file, "ATLAS", atlases[atlas]), header=FALSE)
    n_node <- nrow(network_def)
    n_conn <- n_node * (n_node - 1) / 2

    # initialize empty table for storage of results
    tmp_p_table <- data.frame(Atlas = rep(atlases[atlas], n_test * n_conn),
        Network1 = rep(NA, n_test * n_conn),
        Network2 = rep(NA, n_test * n_conn),
        Node1 = rep(NA, n_test * n_conn),
        Node2 = rep(NA, n_test * n_conn),
        Test_type = rep(test_types, n_conn),
        p = rep(NA, n_test*n_conn))

    # store subject FC for speed
    conn_all <- array(NA, dim=c(n_subj, n_node, n_node))

    for(subj in 1:n_subj) {

        if(!dat$has_run1[subj] & !dat$has_run2[subj]) {
            next
        }

        conn_all[subj, , ] <- conn_matrices[[subj]][[atlas]]

    }

    # mass univariate testing
    conn <- 0
    for(node1 in 1:(n_node - 1)) {

        for(node2 in (node1 + 1):n_node) {

            # Temporary storage of connectivity value
            dat$value <- conn_all[, node1, node2]

            # Fill in table identifier data
            conn <- conn + 1
            tmp_p_table$Node1[((conn - 1) * n_test + 1):(n_test * conn)] <- node1
            tmp_p_table$Node2[((conn - 1) * n_test + 1):(n_test * conn)] <- node2
            tmp_p_table$Network1[((conn - 1) * n_test + 1):(n_test * conn)] <- network_def$V1[node1]
            tmp_p_table$Network2[((conn - 1) * n_test + 1):(n_test * conn)] <- network_def$V1[node2]

            # Remitted greater than HC (t-test, no control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 1] =
                t.test(value ~ group1, dat, alternative="less")$p.value
            # Remitted less than HC (t-test, no control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 2] =
                t.test(value ~ group1, dat, alternative="greater")$p.value

            # Regression to control for covariates
            mod <- lm(value ~ group1 + age + edu + sex + race_bin + site, dat)
            coefs <- coef(summary(mod))
            # Remitted greater than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 3] <- pt(-coefs[2, 3], mod$df)
            # Remitted less than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 4] <- pt(coefs[2, 3], mod$df)

            # any relapse
            mod <- aov(value ~ group3, dat)
            tmp_p_table$p[(conn - 1) * n_test + 5] <- summary(mod)[[1]][["Pr(>F)"]][1]
            mod <- aov(value ~ group3 + age + edu + sex + race_bin + site, dat)
            tmp_p_table$p[(conn - 1) * n_test + 6] <- summary(mod)[[1]][["Pr(>F)"]][1]
            # Pairwise tests
            mod <- dat %>% filter(group3 != "Relapsed") %>%
                lm(value ~ group3 + age + edu + sex + race_bin + site, .)
            coefs <- coef(summary(mod))
            # Remitted greater than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 7] <- pt(-coefs[2, 3], mod$df)
            # Remitted less than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 8] <- pt(coefs[2, 3], mod$df)
            mod <- dat %>% filter(group3 != "Remitted") %>%
                lm(value ~ group3 + age + edu + sex + race_bin + site, .)
            coefs <- coef(summary(mod))
            # Relapsed greater than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 9] <- pt(-coefs[2, 3], mod$df)
            # Relapsed less than HC (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 10] <- pt(coefs[2, 3], mod$df)
            mod <- dat %>% filter(group3 != "Never_Depressed") %>%
                lm(value ~ group3 + age + edu + sex + race_bin + site, .)
            coefs <- coef(summary(mod))
            # Relapsed greater than remitted (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 11] <- pt(coefs[2, 3], mod$df)
            # Relapsed less than remitted (regression, control for covariates)
            tmp_p_table$p[(conn - 1) * n_test + 12] <- pt(-coefs[2, 3], mod$df)

            # survival analysis
            mod <- dat %>%
                filter(group1 == "Remitted") %>%
                coxph(Surv(athf_duration_torelapse, relapsed) ~
                    value + age + edu + sex + race_bin + site, data=.)
            coefs <- coef(summary(mod))
            # Greater FC -> shorter duration to relapse (bad)
            tmp_p_table$p[(conn - 1) * n_test + 13] <- pnorm(coefs[1, 4])
            # Less FC -> shorter duration to relapse (bad)
            tmp_p_table$p[(conn - 1) * n_test + 14] <- pnorm(-coefs[1, 4])

        }

    }

    p_table <- bind_rows(p_table, tmp_p_table)

    # track overall HC tests
    tmp_hc_table <- tibble(Atlas = rep(atlases[atlas], n_test),
        Level = rep("Whole_brain", n_test),
        Network1 = rep(NA, n_test),
        Network2 = rep(NA, n_test),
        Node = rep(NA, n_test),
        Test_type = test_types,
        hc = rep(NA, n_test),
        p = rep(NA, n_test),
        Bonferroni = rep(1, n_test))
    # perform higher criticism
    for(test in 1:n_test) {
        tmp_hc_table$hc[test] <- tmp_p_table %>%
            filter(Test_type == test_types[test]) %>%
            pull(p) %>%
            calc_HC(k1=0.5, emp=TRUE)
    }

    # calculate HC p values (Monte Carlo method)
    tmp_hc_table$p <- calc_HC_p_MC(hc=tmp_hc_table$hc, n_test=n_conn, n_sim=1E5, k1=0.5, emp=TRUE)

    hc_table <- bind_rows(hc_table, tmp_hc_table)

    # perform for network pairs
    for(net1 in 1:n_net) {

        idx1 <- which(network_def$V1 == networks[net1], arr.ind=TRUE)

        for(net2 in 1:net1) {

            cat(paste("   ", networks[net1], "to", networks[net2], "\n"))
            idx2 <- which(network_def$V1 == networks[net2], arr.ind=TRUE)

            # determine number of connections
            if(net1 == net2) {
                n_conn <- length(idx1) * (length(idx1) - 1) / 2
            } else {
                n_conn <- length(idx1) * length(idx2)
            }

            # pull in network-specific results
            tmp_p_table <- p_table %>% filter(Atlas == atlases[atlas],
                ((Network1 == networks[net1] & Network2 == networks[net2]) |
                (Network1 == networks[net2] & Network2 == networks[net1])))

            # track overall HC tests
            tmp_hc_table <- tibble(Atlas = rep(atlases[atlas], n_test),
                Level = rep("Network", n_test),
                Network1 = rep(networks[net1], n_test),
                Network2 = rep(networks[net2], n_test),
                Node = rep(NA, n_test),
                Test_type = test_types,
                hc = rep(NA, n_test),
                p = rep(NA, n_test),
                Bonferroni = rep(NA, n_test))
            # perform higher criticism
            for(test in 1:n_test) {
                tmp_hc_table$hc[test] <- tmp_p_table %>%
                    filter(Test_type == test_types[test]) %>%
                    pull(p) %>%
                    calc_HC(k1=0.5, emp=TRUE)
            }

            # calculate HC p values (Monte Carlo method)
            tmp_hc_table$p <- calc_HC_p_MC(hc=tmp_hc_table$hc, n_test=n_conn, n_sim=1E5, k1=0.5, emp=TRUE)

            hc_table <- bind_rows(hc_table, tmp_hc_table)

            nodes <- sort(unique(c(idx1, idx2)))
            m_node <- length(nodes)

            for(test in 1:n_test) {

                # explore nodes if network is significant
                p_net <- hc_table %>% filter(Atlas == atlases[atlas],
                    Level == "Network",
                    Network1 == networks[net1],
                    Network2 == networks[net2],
                    Test_type == test_types[test]) %>%
                    pull(p)

                # focus on DMN, ECN, and ASN
                # etwork-network FC must be significant after Bonferroni correction
                if(net1 < 4 & net2 < 4 & p_net < (0.05 / 6)) {

                    # track overall HC tests
                    tmp_hc_table <- tibble(Atlas = rep(atlases[atlas], m_node),
                        Level = rep("Node", m_node),
                        Network1 = rep(networks[net1], m_node),
                        Network2 = rep(networks[net2], m_node),
                        Node = nodes,
                        Test_type = rep(test_types[test], m_node),
                        hc = rep(NA, m_node),
                        p = rep(NA, m_node),
                        Bonferroni = rep(m_node, m_node))

                    for(node in 1:m_node) {

                        # pull in node-specific results
                        tmp_p_table <- p_table %>% filter(Atlas == atlases[atlas],
                            ((Network1 == networks[net1] & Network2 == networks[net2]) |
                            (Network1 == networks[net2] & Network2 == networks[net1])),
                            Test_type == test_types[test],
                            (Node1 == nodes[node] | Node2 == nodes[node]))

                        # Remitted greater than HC
                        tmp_hc_table$hc[node] <- calc_HC(p=tmp_p_table$p, k1=0.5, emp=TRUE)
                        n_conn <- nrow(tmp_p_table)

                    }

                    tmp_hc_table$p <- calc_HC_p_MC(hc=tmp_hc_table$hc,
                        n_test=nrow(tmp_p_table), n_sim=1E5, k1=0.5, emp=TRUE)

                    hc_table <- bind_rows(hc_table, tmp_hc_table)

                }

            }

        }

    }

}

