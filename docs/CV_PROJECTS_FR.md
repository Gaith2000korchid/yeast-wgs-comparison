# Descriptions de projets pour les CV ciblés

Ces formulations décrivent les dépôts publics, et non une expérience professionnelle. Sélectionner deux ou trois projets selon l'offre.

## Yeast WGS comparison
**Version courte :** Workflow Snakemake/Python de comparaison WGS Illumina–Nanopore sur levure, avec contrôle qualité, appels de variants et investigation des discordances.

**Version développée :** Construction d'un workflow reproductible sur données publiques appariées : QC, alignement, profondeur, variants et comparaison dans un domaine partagé. Analyse de sensibilité aux filtres ONT et investigation ciblée des discordances avec traçabilité des résultats. Caller exploratoire, sans vérité terrain indépendante.

**Candidatures :** bioinformatique NGS, génomique, ingénierie de recherche, Python scientifique.

## Reproducible cervical biopsy modeling
**Version courte :** Étude prédictive sur données publiques de biopsie cervicale : validation imbriquée groupée, prévention des fuites et analyses de sensibilité.

**Version développée :** Développement d'une analyse Python comparant régression logistique, XGBoost et modèle naïf, avec prétraitement appris dans les plis d'entraînement et quatre sensibilités. Publication des prédictions hors pli, des résultats et d'une reproduction complète ; ROC AUC publiée de 0,613 pour la logistique, concordante au recalcul. Écarts de recalcul XGBoost documentés.

**Candidatures :** data science santé, méthodes quantitatives, biostatistiques prédictives, recherche.

## Biopython M-pro Atlas
**Version courte :** CLI Python/Biopython reliant PROSITE, UniProtKB et PDB/SIFTS pour contrôler la couverture de domaines et les coordonnées structurales.

**Version développée :** Développement d'une chaîne d'interrogation et d'audit d'annotations biologiques : filtrage des structures par couverture du domaine, correspondance des segments SIFTS et analyse des références PROSITE. Sources datées et empreintes, relecture hors ligne et tests Python 3.10/3.12.

**Candidatures :** bioinformatique structurale, développement Python, qualité des données biologiques.

## Viral genome comparator R
**Version courte :** Package R testé et analyse comparative reproductible de 29 génomes de bactériophages, avec alignement MAFFT et analyses taxonomiques.

**Version développée :** Analyse d'un panel public versionné : QC, alignement complet, identité/couverture/distance, conservation et arbre Neighbor-Joining exploratoire. Évaluation quantitative du voisinage taxonomique et des limites du partitionnement en deux clusters ; package et analyses exécutés en CI.

**Candidatures :** génomique comparative, computational biology, développement R, ingénierie de recherche.

## De Bruijn genome assembler
**Version courte :** Assembleur pédagogique Python fondé sur un graphe de De Bruijn pondéré, avec export FASTA, tests de régression et Docker/CI.

**Version développée :** Implémentation des fonctions d'assemblage à partir d'un scaffold universitaire attribué : k-mers, graphes pondérés, simplification des bulles/tips et reconstruction de contigs. Ajout de contrôles FASTQ, de tests ciblés et d'un parcours conteneurisé FASTQ vers FASTA.

**Candidatures :** développement Python, algorithmique, logiciel scientifique, bioinformatique.

## Sélection
- **Vitrine générale :** Yeast WGS, Cervical Biopsy, M-pro Atlas.
- **Génomique/R :** Yeast WGS, Viral Comparator, puis le projet le plus proche de l'offre.
- **Python :** M-pro Atlas, De Bruijn, Yeast WGS.
- **Data science santé :** Cervical Biopsy, puis les expériences biomédicales pertinentes du CV.

Les dépôts ne démontrent pas à eux seuls une plateforme data en production, une validation clinique externe ou un assemblage à grande échelle. Le CV Master reste interne ; ne pas le publier.
