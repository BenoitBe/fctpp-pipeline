# Audit du code R — FC.TPP et `rDEP_Version_Definitive_20625.R`

**Date :** 27/09/2026 — la version Python, d'abord appelée `fctpp`, a été renommée **Olympus** le 17/11/2026.

**Périmètre :**
- package R `FC.TPP` 0.0.0.9000 (archive `FC.TPP.zip`, index git sans commit) ;
- script `rDEP_Version_Definitive_20625.R` (V1.0 du 20/06/2025).

**Méthode :**
1. Lecture du code.
2. Réimplémentation Python, puis comparaison avec les sorties R livrées (`TPP_results.xlsx`).
3. Recalcul depuis les `T*.tsv` bruts.
4. Exécution réelle de Deimos (commit `4d49d24`) sur les 12 températures, refaite le 17/11/2026 avec `e160463` : mêmes résultats.

Pas de R disponible dans l'environnement : les constats sur DEP sont vérifiés sur les sorties R existantes.

---

## 1. Synthèse

| # | Sévérité | Constat | Dans la version Python (Olympus) |
|---|---|---|---|
| A1 | Élevée | Appariement échantillons ↔ design par effet de bord (`delete_prefix` de DEP) | Fait en amont par Deimos : un design par température, appariement exact des labels |
| A2 | Élevée | « 2 températures consécutives » = 2 *lignes successives* : trous et inversions de signe acceptés | Règle `strict` par défaut ; `r_legacy` pour reproduire R |
| A3 | Élevée | `p.adj` de DEP = q-value fdrtool, **inférieure à p.val dans 29 à 47 % des lignes** | Deimos : BH ; Olympus utilise `p.val` (comme R) ou `p.adj` BH au choix |
| A4 | Élevée | Aucun contrôle de multiplicité entre températures × contrastes | Documenté ; règle stricte ; pistes de calibration |
| A5 | Élevée | Imputation MinProb non seedée et déterminante pour certains appels | Drapeau `imputed_in_run` ; imputation déléguée à Deimos (QRILC, seed) |
| A6 | Moyenne à élevée | Abondance et stabilité confondues (aucune référence basse température) | Drapeau `shift_class = abundance-like` |
| B1 | Moyenne | Les exports Excel s'écrasent d'une fonction à l'autre ; échec au 2e run | Un classeur unique, écrit une fois |
| B2 | Moyenne | Températures codées en dur `T1..T12` | Série lue dans les données ; °C optionnels |
| B3 | Moyenne | Trois jeux de seuils incohérents ; la doc annonce p ajustée, le code utilise p brute | Seuils volcano = seuils d'appel par défaut ; doc alignée |
| B4 | Moyenne | Dossier par défaut `tpp_data` au lieu de `TPP_data` (Linux/HPC) | Recherche insensible à la casse |
| B5 | Faible | Message « proteins were imputed » qui compte en réalité les protéines filtrées | Sans objet (filtrage dans Deimos, Methods explicite) |
| B6 | Moyenne | `.Rbuildignore` exclut les données d'exemple du package construit | Sans objet (données réelles non versionnées) |
| B7 | Faible | 3D : 1 trace par protéine (5,6 Mo/HTML), plans en dur, PDF non fermé sur erreur | 1 trace de fond, plans calculés, gestionnaires de contexte |
| C1 | Faible | Hygiène du dépôt (sorties versionnées, doublons de données, README vide, test factice) | Nouveau dépôt, CI, 70 tests |

---

## 2. Constats détaillés — package FC.TPP

### A1. Appariement échantillons ↔ design implicite

`find_design()` lit un seul design dont les labels sont ceux de **T1** (`Sple.1B` … `Sple.1J`). Pour T2 à T12, les colonnes s'appellent `Sple.2B` … `Sple.12J`. `DEP::make_se()` apparie quand même, parce qu'il retire le plus long préfixe commun des deux côtés avant de comparer (`delete_prefix`). Seule la lettre finale (B à J) est donc réellement comparée.

Deux dépendances implicites en découlent :
- à la conversion `-` → `.` de `read.delim(check.names = TRUE)`, dont dépend aussi la regex `Sple\\.` de `format_tsv` ;
- au fait que les labels d'une température ne partagent pas davantage de préfixe que ceux de T1.

Rien n'est rapporté à l'utilisateur.

*Vérification.* La même logique, réimplémentée explicitement, redonne exactement les présences R. Sur les 12 températures, T1 s'apparie en exact et T2 à T12 par suffixe.

*Dans la nouvelle chaîne.* L'appariement se fait dans Deimos, qui exige un design par température avec les labels exacts de cette température (`Sple-2B` à T2, etc.) et arrête le run si un label manque. Le problème disparaît donc, à condition de fournir X designs cohérents : même ordre des conditions et même contrôle partout. Olympus vérifie ensuite que les X classeurs exposent les mêmes contrastes.

### A2. Règle des températures consécutives (`tidy_tpp`, l. 989-1000)

```r
cond_k      <- abs(diff) >= FC_cutoff & logpval >= pval_cutoff
cond_suiv_k <- c(cond_k[-1], FALSE)
seuil_k     <- any(cond_k & cond_suiv_k)
```

Les lignes sont celles où la protéine est présente, triées par température. Deux conséquences :
- **les trous sont ignorés** : absente à T5, la protéine voit T4 et T6 traitées comme consécutives ;
- **le signe n'est pas vérifié** : +0,8 à T4 puis −0,8 à T5 donne « Over ».

*Impact mesuré sur `TPP_results.xlsx`.*
- Aux seuils par défaut (FC 0,6 ; p < 0,01), aucun impact : 2 hits, tous deux adjacents et de même signe.
- À FC 0,3 et p < 0,05, **6 des 29 appels** reposent uniquement sur une paire non adjacente ou de signe inversé :

  | Contraste | Appels concernés |
  |---|---|
  | C1_vs_Ctl | 3/18 |
  | C2_vs_C1 | 1/7 |
  | C2_vs_Ctl | 2/4 |

*Python.*
- `rule: strict` (défaut) : températures adjacentes dans la série, même signe, `min_consecutive` paramétrable.
- `rule: r_legacy` : reproduit R à l'identique (validé sur 23 792 lignes).
- Les deux verdicts figurent dans la feuille `Calls` (`over_strict`, `over_r_legacy`).

### A3. `p.adj` de DEP n'est pas une correction BH

`DEP::test_diff()` calcule `p.adj` avec `fdrtool` sur les statistiques t. fdrtool estime une distribution nulle *empirique* ; ce n'est pas une correction de p-values. Dans le `TPP_data` livré :

| Contraste | p.adj < p.val | p.adj < 0,05 | p.val < 0,05 |
|---|---|---|---|
| C1_vs_Ctl | 30,9 % des lignes | 599 | 287 |
| C2_vs_C1 | 47,3 % | 901 | 289 |
| C2_vs_Ctl | 29,4 % | 559 | 269 |

À T1 (C1_vs_Ctl), on compte **96** protéines à p.adj < 0,05 contre **10** à p.val < 0,05. Exemple : Q2FVX1 à T1 a p = 0,137 et « p.adj » = 0,006.

Conséquences :
- **FC.TPP.** Les appels utilisent `p.val`, donc les résultats TPP ne sont pas touchés. En revanche, les docstrings de `magma()` et de `tidy_tpp()` annoncent une p-value ajustée.
- **rDEP.** `add_rejections(alpha = 0.05, lfc = log2(1.2))` construit la colonne `significant` sur ces q-values. **Les listes « significatives » produites par rDEP sont probablement très optimistes.** Il faut les revoir si elles ont servi à des interprétations.

*Python.* Deimos fournit une correction BH par contraste. Olympus garde `p.val` par défaut pour la parité, et `call.p_column: p.adj` permet d'appeler sur BH.

### A4. Multiplicité

On compte environ 2 041 protéines × 12 températures × 3 contrastes ≈ 73 000 tests, chacun seuillé à p brut < 0,01.

La contrainte de 2 températures consécutives réduit les faux positifs, mais sans les calibrer. Les températures adjacentes proviennent en effet des mêmes lysats (même réplicat biologique) : elles sont positivement corrélées, et le taux de faux positifs est donc supérieur à ce que suggérerait l'indépendance.

*Pistes.*
- Permutations des labels de condition au sein des réplicats, pour un FDR empirique de la règle complète.
- Approche par courbes, de type NPARC (Childs *et al.*, 2019).
- Modèle limma unique avec interaction température × condition.

### A5. Imputation

`impute(fun = "MinProb", q = 0.01)` est appelé sans `set.seed` : deux exécutions ne donnent pas les mêmes p-values pour les protéines imputées. C'est le même constat que dans l'audit plateforme du 22/09.

Surtout, à haute température, l'absence d'une protéine dans une condition est souvent **vraie** (précipitation : MNAR). La valeur imputée fixe alors l'amplitude du FC.

*Mesures.*
- Sur les lignes sans imputation (22 423), les log2FC Deimos et DEP sont **identiques à 2 × 10⁻¹⁴**.
- Sur les lignes imputées (1 369), l'écart médian est de 0,29 log2 (QRILC contre MinProb).
- Le seul hit supplémentaire de la version Deimos (Q2FVV0, C2_vs_Ctl, T2–T3) repose sur une température imputée et sur une protéine absente à T1.

*Python.* Drapeau `imputed_in_run` sur chaque hit, à partir de `num_NAs` de Deimos. Ce compte porte sur la protéine toutes conditions confondues ; le drapeau est donc conservateur. Imputation et seed sont gérés par Deimos, et leur homogénéité entre températures est vérifiée.

### A6. Abondance et stabilité thermique

L'analyse compare les conditions **à chaque température séparément**. Un FC observé dès la température la plus basse traduit une différence d'abondance (expression, solubilité), pas un déplacement de Tm. Or une protéine simplement surexprimée dans C1 franchit le seuil à toutes les températures, et R la déclare hit.

*Python.* `shift_class = abundance-like` quand le hit passe déjà les seuils, avec le même signe, à la température de référence (la plus basse). C'est une heuristique : la solution propre est le modèle d'interaction évoqué en A4.

Sur les données, les 2 hits R sont `thermal` (T9–T10, FC ≈ +1 à +1,8 ; transporteurs biotine et riboflavine, C1_vs_Ctl).

### B1. Exports Excel qui s'écrasent

`format_tsv`, `sumXP`, `process_and_export` et `tidy_tpp` écrivent tous par défaut dans `TPP_results.xlsx`. `tidy_tpp()` recrée le classeur (`createWorkbook` + `overwrite = TRUE`) : les feuilles `comparaison_T*` et `-log2_*` sont perdues. **C'est confirmé** : le `TPP_results.xlsx` livré ne contient que `TPP_data` et `Stats_proteins`.

À l'inverse, `process_and_export(overwrite_excel = FALSE)` recharge le fichier existant puis appelle `addWorksheet` sur un nom déjà présent : erreur dès le 2e run.

### B2. Série de températures codée en dur

- `tidy_tpp` : `factor(..., levels = paste0("T", 1:12))`.
- `plot_and_3D` : plans de seuil sur `x = c("T1", "T12")`.

Un autre nommage (T0, 37C…) ou un autre nombre de températures donne des `NA` silencieux, puis des protéines écartées. Enfin, aucune valeur en °C n'est utilisée (axe catégoriel).

### B3. Seuils incohérents

| Étape | FC | p |
|---|---|---|
| `add_rejections` | log2(1,2) = 0,26 | p.adj (fdrtool) < 0,05 |
| `magma` (volcano) | 0,3 | p.val < 0,01 |
| `tidy_tpp` / `plot_and_3D` | 0,6 | p.val < 0,01 |

Les volcanos ne montrent donc pas ce qui est appelé.

### B4 à B7

- **B4.** `find_design(folder = "tpp_data")` alors que les données sont dans `TPP_data` : échec sur système de fichiers sensible à la casse (Linux, cluster).
- **B5.** `process_and_export` affiche « *n* proteins were imputed » avec *n* = protéines avant filtrage − protéines après. Ce sont des protéines **retirées**, pas imputées.
- **B6.** `.Rbuildignore` : le motif non ancré `TPP_data/` et le motif `^.*\.xlsx$` excluent `inst/extdata/TPP_data/*` (y compris `experimental_design.xlsx`) du tarball. Le package construit n'a plus ses données d'exemple.
- **B7.** `plot_and_3D` :
  - une trace plotly par protéine (`split = ~name`), d'où 5,6 Mo par HTML ;
  - nom de fichier construit avec la dernière valeur de `pick` (« All_prot » même quand le graphe contient tout) ;
  - `pdf()` ouvert sans `on.exit(dev.off())` ;
  - `size =` déprécié depuis ggplot2 3.4 (`linewidth`) ;
  - vue 2D « All_prot » : ~1 500 couleurs sans légende.

  Autres points mineurs :
  - `sumXP` nomme ses feuilles `-log2_` alors que les valeurs sont en log2, et son message d'erreur colle le data.frame entier ;
  - `format_tsv` ne retire les contaminants que si le groupe *commence* par `cRAP` (groupes mixtes conservés) ;
  - `ID` et `name` sont inversés par rapport à la convention DEP : les figures sont étiquetées par accession plutôt que par nom de gène.

### C1. Hygiène du dépôt

- Index git sans aucun commit, mais **~32 Mo de sorties indexées** : 3 HTML de 5,6 Mo, `TPP_plots.pdf` 7 Mo, `volcano_plots.pdf` 4 Mo, `TPP_results.xlsx` 4 Mo.
- Données en double (`TPP_data/` et `inst/extdata/TPP_data/`).
- `.RData` de 12,7 Mo, `.Rhistory` et `.Rproj.user` dans l'archive transmise.
- README issu du template, non rempli ; seul test : `expect_equal(2 * 2, 4)`.
- `magick` et `webshot2` (Chrome) en `Imports` pour une option GIF.

---

## 3. Script `rDEP_Version_Definitive_20625.R`

C'est l'ancêtre mono-jeu de `process_and_export`. **Deimos en est l'équivalent Python maintenu** : même regex `Sple`, même règle cRAP, `filter_missval(thr = 1)`, limma `~0 + condition`, ANOVA et heatmap. Il n'a donc pas été reporté une seconde fois. Constats :

| Point | Détail | Couvert par Deimos |
|---|---|---|
| Significativité | `add_rejections` sur les q-values fdrtool (cf. A3) : listes trop permissives | Oui (BH) |
| Normalisation | `normalize_vsn()` calculé mais **jamais utilisé** (`data_norm` orphelin) ; `plot_normalization` n'affiche que les données filtrées | Oui (pas de renormalisation, diagnostic RLE) |
| Design | Codé en dur (12 échantillons WT/Mutant) ; contrôles de doublons exécutés *après* `make_se` ; `tsv_brut$Protein.IDs` n'existe pas | Oui (CSV + diagnostic des labels) |
| Imputation | MinProb sans seed | Oui (QRILC/auto + seed, correctifs du 22/09) |
| Volcanos | Seuils en dur (p brute 0,01 ; FC 0,3) | Oui (YAML) |
| ANOVA | Boucle `rbind` quadratique ; `aov` sur valeurs imputées ; BH correct | Oui |
| Heatmap | `gpar(fontsize = 0,25)` : virgule décimale, donc police 0 (noms de lignes invisibles) et `25` passé sans nom | Oui |
| Excel | Feuille `Plots` vide (`insertImage` commentés, contrastes d'un autre projet) | Oui |
| GO (topGO) | Voir le détail ci-dessous | Partiellement (gProfiler / STRING / annotation locale) |

Problèmes de l'enrichissement GO (topGO) :
- GAF d'un autre organisme codé en dur (*S. agalactiae*) ;
- sélection sur p brute ≤ 0,01, sans critère de FC ni de sens (up et down mélangés) ;
- gènes nommés par `Genes` (vides ou multiples) ;
- p-values `NA` propagées dans `allGenes` ;
- `GenTable` renvoie « < 1e-30 » en texte, que `as.numeric` transforme en `NA` ;
- aucune correction multiple.

---

## 4. Validation de la version Python (Olympus)

| Contrôle | Résultat |
|---|---|
| Appels (`r_legacy`) sur `TPP_data` R | 23 792 lignes, 3 contrastes : `_logpval` et `_threshold` **identiques** |
| Statistiques de présence recalculées depuis les `T*.tsv` bruts | 2 215 protéines : `nb_temp` et `max_consecutive` **identiques** ; 2 041 retenues (idem R) |
| 12 runs Deimos (un par température) puis `olympus analyze` sur les 12 classeurs | 12/12 runs OK (~30 s chacun) ; mêmes 2 041 protéines retenues |
| log2FC Deimos vs DEP, lignes sans imputation | identiques (écart max 2 × 10⁻¹⁴) |
| Hits | Les 2 hits R retrouvés + 1 hit supplémentaire signalé comme douteux (`imputed_in_run`, absent à T1) |
| Tests | 70 tests (synthétiques + CLI) ; `-m realdata` et `-m deimos` sur les données réelles |

---

## 5. Recommandations hors FC.TPP

1. **Deimos, `load_data`** : `tsv.columns[:6]` suppose 6 colonnes d'annotation. Un export DIA-NN ancien (4 colonnes) fait basculer 2 runs en annotations ; le diagnostic des labels s'arrête alors avec un message trompeur. Détecter plutôt les colonnes numériques. En attendant, exporter les matrices TPP avec DIA-NN ≥ 2 (6 colonnes d'annotation).
2. **Deimos, regex `Sple-[^_]+`** codée en dur (diagnostic et chargement) : à passer en paramètre YAML.
3. **Deimos, seed** : les correctifs du 22/09 (seed, QRILC par régression quantile, `impute_method: auto`) ne sont pas dans le commit GitHub `4d49d24` utilisé ici. *Mise à jour du 17/11 :* le seed est dans `main` (`e160463`, `seed: 42` par défaut). Restent ouverts les points 1 et 2.
4. **Deimos en mode TPP** : éviter `impute_method: auto`, qui peut choisir une méthode différente selon la température. Fixer la même méthode pour les X runs ; `olympus analyze` vérifie l'homogénéité dans les feuilles `Methods` et la signale sinon.
5. **rDEP** : revoir toute liste de protéines « significatives » produite avec `add_rejections` (cf. A3).
