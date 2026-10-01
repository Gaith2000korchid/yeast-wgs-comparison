# Synthèse du projet : comparaison WGS Illumina / Nanopore chez la levure

## Question et démarche

Comment construire une comparaison reproductible des SNPs obtenus avec deux technologies de séquençage, et expliquer les limites des résultats ? Ce projet traite les lectures appariées Illumina et les longues lectures Oxford Nanopore d'un même isolat de *Saccharomyces cerevisiae*, CIC/Ponton11 (P11), échantillon ENA SAMEA3895683. La référence est S288C R64-1-1, Ensembl 113, avec un modèle nucléaire diploïde documenté.

Le workflow **Snakemake** relie les fichiers d'entrée aux rapports : FastQC/fastp pour Illumina, NanoPlot pour ONT, alignement BWA-MEM/minimap2, SAMtools pour les BAMs et la couverture, BCFtools pour l'appel, la normalisation et le filtrage, puis comparaison des SNPs alternatifs dans un masque de profondeur partagé. MultiQC rassemble les rapports de QC.

Les données réelles comprennent **2 105 894 paires Illumina**, sélectionnées de manière déterministe à 25 %, et **29 914 lectures ONT consensus 2D** historiques, pass et fail. Les records template/complement corrélés d'une même molécule ne sont pas cumulés. Les versions, paramètres, checksums, logs et entrées sont documentés dans le [protocole de reproduction](P11_RUN.md).

## Résultats du génome entier

| Indicateur | Illumina | ONT historique |
|---|---:|---:|
| Lectures primaires alignées | 97,40 % | 74,64 % |
| Profondeur nucléaire moyenne filtrée, BQ≥13 | 45,93× | 1,72× |
| Profondeur nucléaire moyenne, ONT BQ≥7 | identique | 8,43× |

La qualité moyenne ONT est Q8, avec un N50 de 7 754 bases. Le seuil BQ13 retient donc peu de couverture utilisable. Les alertes FastQC de contenu GC sont conservées et décrites ; aucune analyse taxonomique de contamination n'a été réalisée.

| Comparaison nucléaire | ONT BQ13 | ONT BQ7 |
|---|---:|---:|
| Bases dans le masque de profondeur partagé | 96 235 (0,80 %) | 6 435 340 (53,31 %) |
| SNPs alternatifs partagés | 156 | 5 444 |
| SNPs Illumina-only | 140 | 32 968 |
| SNPs ONT-only | 4 | 27 |
| Jaccard des SNPs alternatifs | 52,00 % | 14,16 % |
| Accord des génotypes aux seuls SNPs partagés | 98,08 % | 98,68 % |

Seul le seuil de qualité des bases ONT change entre ces deux analyses. Les deux Jaccard portent sur des domaines différents : leur différence ne mesure pas une amélioration ou dégradation de précision. L'accord élevé des génotypes partagés est fortement conditionné par les sites retenus. À BQ7, le caller ONT n'émet que huit SNPs nucléaires hétérozygotes parmi 5 545 records alternatifs. Ce baseline sous-représente donc les appels hétérozygotes.

Les [résultats complets et leurs preuves](results/P11/README.md) permettent de vérifier chaque nombre.

## Revue de lectures et audit régional

Quinze loci ont été sélectionnés de manière déterministe pour mesurer leur support diagnostique ; **six** ont ensuite été examinés par Gaith Korchid dans IGV Web. Les [deux captures archivées et les comptes transmis](results/P11/review/IGV_WEB_REVIEW_FR.md) montrent notamment :

- **III:105658 A>T** : Illumina affiche 59 A et aucun T ; ONT affiche 7 T et 3 A. Discordance non résolue.
- **IV:308249 T>C** : Illumina affiche 69 T ; ONT affiche 12 C et 1 A, avec cinq délétions signalées. Le contexte ONT comporte des gaps et insertions à examiner.
- **I:2533 G>T** : mélange G/T en Illumina ; uniquement T parmi huit bases ONT observées. Cela ne valide pas une homozygotie ONT.
- **I:27486, IV:1364942 et IX:37309** : l'allèle alternatif est soutenu dans les deux plateformes malgré le classement ONT-only des VCF filtrés.

Pour ces trois derniers sites, une [reconstruction régionale du caller](results/P11/review/REGIONAL_CALLER_AUDIT_FR.md) appelle les SNPs Illumina avant de les retirer au filtrage. FORMAT/DP vaut 3, 2 et 6, sous le seuil 8 ; IV:1364942 échoue aussi à QUAL≥30. Désactiver uniquement BAQ augmente FORMAT/DP de 2 à 26 et de 6 à 46 aux deux derniers sites, qui passent alors les filtres. Le premier reste à profondeur 3.

Cette expérience distingue **support dans les lectures**, **appel d'un variant** et **conservation après filtrage**. Ses valeurs proviennent des BAMs régionaux ; elles ne remplacent pas les intermédiaires originaux du run génome entier. BAQ désactivé est un contrôle diagnostique et n'est pas adopté comme nouveau réglage final.

## Ce qui est démontré et ce qui reste ouvert

Le projet démontre l'exécution reproductible d'un workflow à deux branches, la traçabilité des données, le contrôle des coordonnées et dénominateurs, une expérience de sensibilité et une investigation ciblée d'une discordance de VCF. **21 tests Python** et un workflow synthétique complet vérifient le logiciel, y compris la récupération de trois SNPs plantés, les données appariées, les masques et l'export régional.

Il n'établit pas une précision biologique, une supériorité de plateforme ni une concordance des génotypes à l'échelle du génome. Il manque une vérité indépendante, des exclusions validées de répétitions/mappabilité et un caller adapté à la chimie ONT historique. Les neuf autres loci sélectionnés ne sont pas examinés visuellement. Les réglages IGV exacts sont inconnus. La préparation ONT consensus 2D ne représente pas les performances des données ONT modernes. Les indels et variants structuraux ne sont pas comparés.

## Présentation orale en une minute

« J'ai construit un workflow Snakemake comparant les données Illumina et Nanopore d'un même isolat de levure. J'ai traité les lectures, aligné les deux technologies et comparé les SNPs dans les régions couvertes par les deux jeux de données. Les anciennes lectures Nanopore avaient une qualité moyenne Q8 : le seuil BQ13 limitait la comparaison à 0,8 % du génome nucléaire. Une expérience BQ7 a élargi ce domaine à 53,3 %, sans permettre d'affirmer une précision biologique. J'ai ensuite examiné six sites dans IGV. Trois variants classés ONT-only avaient aussi du support Illumina. Un audit régional a montré leur retrait par les filtres du caller et une forte sensibilité à BAQ pour deux sites. Le résultat principal est donc une comparaison traçable, avec une explication concrète de l'effet des paramètres et des limites de l'interprétation. »

Implémentation développée avec assistance IA ; l'examen IGV et les comptes transmis sont attribués à Gaith Korchid. Les sources et outils tiers sont crédités dans le README.
