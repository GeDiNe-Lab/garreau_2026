library(Matrix)
library(tidyverse)
library(parallel)
library(RANKS)

# ============================================================
# build_permutation_bank()
# ------------------------------------------------------------
# Construit UNE FOIS une banque de profils RWR sur seeds aleatoires.
# A reutiliser ensuite pour normaliser n'importe quel profil unique
# ou groupe de patients, sans jamais relancer RWR.
#
# @param RWR_net        matrice d'adjacence du reseau (genes en
#                        rownames/colnames)
# @param seed_sizes      vecteur des tailles de seed-sets a tirer
#                        (ex: le nombre de seeds observe chez tes
#                        vrais patients -> reflete la vraie distribution).
#                        Un seul entier = taille fixe pour tous les tirages.
# @param n_permutations  nombre de profils aleatoires a generer (taille
#                        de la banque). Parametre unique controlant le
#                        nombre de permutations calculees.
# @param gamma, tmax, eps, norm  parametres RWR (identiques a ta fonction RWR())
# @param mc.cores        nombre de coeurs (mclapply si > 1)
# @param verbose         afficher la progression
#
# @return matrice numerique genes x n_permutations, avec un attribut
#         "seed_sizes" (taille du seed-set utilise pour chaque colonne)
# ============================================================
build_permutation_bank <- function(RWR_net, seed_sizes,
                                   n_permutations = 1000,
                                   gamma = 0.7, tmax = 1000, eps = 1e-10, norm = TRUE,
                                   mc.cores = 1, verbose = TRUE) {
    all_genes <- rownames(RWR_net)

    run_one <- function(i) {
        size_i <- sample(seed_sizes, 1)
        random_seeds <- sample(all_genes, size_i)

        if (verbose && i %% 100 == 0) message("Permutation ", i, "/", n_permutations)

        score <- RWR(RWR_net, random_seeds,
            gamma = gamma, tmax = tmax, eps = eps, norm = norm
        )$p %>%
            unname()

        list(score = score, size = size_i)
    }

    if (mc.cores > 1) {
        results <- mclapply(seq_len(n_permutations), run_one, mc.cores = mc.cores)
    } else {
        results <- lapply(seq_len(n_permutations), run_one)
    }

    bank <- sapply(results, function(x) x$score)
    rownames(bank) <- all_genes
    colnames(bank) <- paste0("perm_", seq_len(n_permutations))

    attr(bank, "seed_sizes") <- sapply(results, function(x) x$size)

    return(bank)
}


# ============================================================
# normalize_rwr_score()
# ------------------------------------------------------------
# Normalise (Z-score) un profil RWR unique OU une matrice de profils
# (auquel cas la moyenne des patients est calculee puis normalisee).
#
# Le null est construit en reechantillonnant/agregant la banque de
# permutations precalculee, a la taille N du profil fourni (N = 1
# pour un profil unique, N = ncol(score_input) pour une matrice).
# -> aucune propagation RWR n'est relancee si permutation_bank est fourni.
#
# @param score_input      vecteur (profil unique) OU matrice genes x
#                         patients (plusieurs profils a moyenner)
# @param RWR_net          reseau associe (requis seulement si
#                         permutation_bank est NULL, pour construire
#                         la banque a la volee)
# @param permutation_bank banque precalculee (sortie de
#                         build_permutation_bank). Si NULL, elle est
#                         construite ici (recalcul complet).
# @param n_permutations   si permutation_bank est NULL : taille de la
#                         banque a construire. Si permutation_bank est
#                         fourni : nombre de tirages nuls (bootstrap)
#                         utilises pour estimer moyenne/sd du null a
#                         la taille N. Parametre unique pour fixer le
#                         nombre de permutations dans les deux cas.
# @param seed_sizes       requis seulement si permutation_bank est NULL
#                         (voir build_permutation_bank)
# @param gamma, tmax, eps, norm, mc.cores  transmis a build_permutation_bank
#                         si la banque doit etre construite
# @param verbose          afficher la progression
#
# @return vecteur nomme (genes) de scores Z-normalises
# ============================================================
normalize_rwr_score <- function(score_input, RWR_net = NULL,
                                permutation_bank = NULL,
                                n_permutations = 1000,
                                seed_sizes = NULL,
                                gamma = 0.7, tmax = 1000, eps = 1e-10, norm = TRUE,
                                mc.cores = 1, verbose = TRUE) {
    # --- 1. Profil unique vs matrice : calcule le profil a normaliser et N ---
    if (is.null(dim(score_input))) {
        profile <- score_input
        N <- 1
    } else if (ncol(score_input) == 1) {
        profile <- score_input[, 1]
        N <- 1
    } else {
        profile <- rowMeans(score_input, na.rm = TRUE)
        N <- ncol(score_input)
    }

    # --- 2. Recupere ou construit la banque de permutations ---
    if (is.null(permutation_bank)) {
        if (is.null(RWR_net) || is.null(seed_sizes)) {
            stop(
                "RWR_net et seed_sizes sont requis pour construire la banque ",
                "si permutation_bank n'est pas fourni"
            )
        }
        if (verbose) {
            message(
                "Aucune banque fournie : construction de ", n_permutations,
                " permutations (recalcul complet)..."
            )
        }
        permutation_bank <- build_permutation_bank(
            RWR_net,
            seed_sizes = seed_sizes, n_permutations = n_permutations,
            gamma = gamma, tmax = tmax, eps = eps, norm = norm,
            mc.cores = mc.cores, verbose = verbose
        )
        n_null_draws <- n_permutations
    } else {
        if (verbose) {
            message(
                "Banque de permutations reutilisee (taille ",
                ncol(permutation_bank), "), pas de recalcul RWR."
            )
        }
        n_null_draws <- n_permutations
    }

    B <- ncol(permutation_bank)
    replace_flag <- N > B
    if (replace_flag && verbose) {
        message("N (", N, ") > taille de la banque (", B, ") : tirage avec remise.")
    }

    # --- 3. Reechantillonne/agrege la banque a la taille N pour batir le null ---
    null_means <- replicate(n_null_draws, {
        idx <- sample(seq_len(B), size = N, replace = replace_flag)
        if (N == 1) permutation_bank[, idx] else rowMeans(permutation_bank[, idx, drop = FALSE])
    })
    # null_means : genes x n_null_draws

    null_mean_per_gene <- rowMeans(null_means)
    null_sd_per_gene <- apply(null_means, 1, sd)

    z_profile <- (profile - null_mean_per_gene) / null_sd_per_gene
    names(z_profile) <- rownames(RWR_net) %||% rownames(permutation_bank)

    return(z_profile)
}


# ============================================================
# Exemple d'utilisation
# ============================================================
# # 1. Construire la banque UNE FOIS (ex: 2000 permutations, tailles de
# #    seeds tirees de la distribution reelle observee chez tes patients)
# observed_seed_counts <- sapply(seeds_list, length)
# bank <- build_permutation_bank(
#     RWR_net, seed_sizes = observed_seed_counts,
#     n_permutations = 2000, mc.cores = 4
# )
# save(bank, file = "article/transit_data/permutation_bank.RData")
#
# # 2. Reutiliser la banque pour normaliser un patient unique
# z_patient <- normalize_rwr_score(RWR_results[, "Patient_1"], permutation_bank = bank)
#
# # 3. Reutiliser la MEME banque pour normaliser un groupe de 20 patients
# z_group20 <- normalize_rwr_score(RWR_results[, sample(colnames(RWR_results), 20)],
#                                   permutation_bank = bank)
