# Démarrer et comprendre le projet

## Première étape : exécuter le petit test

Sous Windows, ouvre Ubuntu avec WSL2. Suis les commandes du README pour installer l'environnement et lancer le test synthétique. Le générateur crée une petite séquence aléatoire, trois SNP connus et deux types de lectures. Ce n'est pas un génome de levure et ce test ne mesure pas la performance biologique du pipeline.

Ouvre le rapport HTML. Vérifie où sont stockés les BAM, les VCF, les métriques et les logs. Relance le workflow : Snakemake doit éviter de recalculer les sorties déjà à jour.

## Comprendre une branche à la fois

1. **FastQC / fastp** : observer la qualité des lectures Illumina, retirer les adaptateurs et filtrer/nettoyer les lectures.
2. **BWA / minimap2** : retrouver la position des lectures sur la référence ; les technologies utilisent des aligneurs différents.
3. **SAMtools** : trier/indexer les alignements, marquer les doublons Illumina et mesurer la couverture.
4. **BCFtools** : chercher des différences par rapport à la référence, normaliser leur représentation et appliquer les filtres.
5. **Comparaison** : compter les SNP communs uniquement dans les zones où les deux jeux de lectures ont une couverture acceptable.

Pour Nanopore, bcftools sert ici de première base exploratoire. Nous devrons justifier un caller adapté aux données réelles et documenter ses limites avant d'interpréter les différences comme biologiques.

## Passage aux données réelles

Le document DATASETS.md identifie une étude publique de levure avec Illumina et ONT. La prochaine étape scientifique consiste à sélectionner un seul isolat, confirmer que les deux technologies correspondent au même échantillon, vérifier la ploïdie, récupérer les FASTQ et leur provenance, puis exécuter le workflow.

Le premier livrable scientifique sera un tableau de QC, les statistiques de SNP, et quelques discordances examinées dans IGV. Une bonne explication doit dire ce qui a été observé, comment le pipeline le calcule et ce qui reste incertain.

## Ce que tu dois pouvoir expliquer en entretien

- Pourquoi utiliser BWA pour Illumina et minimap2 pour Nanopore ?
- Pourquoi la ploïdie change-t-elle l'appel des génotypes ?
- Pourquoi comparer les mêmes régions couvertes par les deux technologies ?
- Pourquoi un SNP présent dans un seul VCF n'est-il pas automatiquement faux ?
- Comment reproduire une exécution et identifier l'étape qui a échoué ?
- Qu'as-tu développé personnellement à partir des tutoriels et sources cités ?

Le projet est développé avec assistance IA. Lis chaque étape, exécute-la et vérifie les résultats pour pouvoir défendre réellement le travail.
