# Inspection visuelle P11 dans IGV : guide pratique

Le but est de décrire le support des lectures, les différences d'alignement et le contexte local. **Aucune inspection manuelle dans IGV n'est encore déclarée terminée.** Une observation visuelle ne constitue pas à elle seule une validation biologique.

## Récupérer le paquet

Ouvrir l'onglet [Actions du dépôt](https://github.com/Gaith2000korchid/yeast-wgs-comparison/actions/workflows/igv-bundle.yml), puis une exécution réussie de **Export real P11 IGV bundle**. En bas de la page, télécharger l'artifact **P11-IGV-review**. GitHub demande une connexion pour télécharger les artifacts. Leur conservation est de 90 jours ; ensuite utiliser **Run workflow** pour reconstruire le paquet.

Sur Windows, extraire l'artifact, puis extraire `igv_bundle.zip` dans un dossier, par exemple `Documents\yeast-igv`. Le dossier `igv_bundle` doit contenir ensemble `session.xml`, `reference.fa`, `reference.fa.fai`, les deux `.bam` et leurs `.bai`. Ne pas les ouvrir directement dans l'archive ZIP.

Le paquet contient la référence complète et les alignements chevauchant des fenêtres de ±2 000 bases autour des 15 loci sélectionnés. Les enregistrements originaux sont conservés entiers, sans filtrage de qualité lors de l'extraction et sans dupliquer un alignement traversant plusieurs fenêtres. Le script vérifie que les comptages aux 15 loci sont identiques entre les BAMs complets et régionaux. Le workflow vérifie aussi l'identité avec les mesures déjà publiées.

**Hors des fenêtres, la couverture est incomplète**, même si des lectures peuvent être visibles. Les mates situés entièrement hors des fenêtres peuvent être absents. Ne pas utiliser ce paquet pour la couverture du génome entier, un nouvel appel de variants ou un bilan global des paires. Le fichier `extraction_windows.bed` donne le périmètre exact.

## Ouvrir IGV

1. Télécharger [IGV Desktop pour Windows avec Java inclus](https://igv.org/doc/desktop/DownloadPage/) et l'installer.
2. Dans IGV, choisir **File → Open Session** et sélectionner `igv_bundle/session.xml`. La référence et les deux pistes BAM doivent apparaître. Le XML et ses chemins sont testés ; l'ouverture dans l'interface graphique reste à vérifier sur le PC du reviewer.
3. Saisir d'abord `III:105558-105758` dans la zone de navigation. La position du SNP est **105658**. Garder ce nom de chromosome romain : `III`, pas `chr3`.

## Les deux loci prioritaires

Les fractions ci-dessous sont des comptages externes à IGV, avec MAPQ≥20, BQ≥13 pour Illumina et BQ≥7 pour ONT. Les mates Illumina qui se chevauchent comptent une seule fois. IGV peut afficher une profondeur différente selon ses filtres et son sous-échantillonnage.

| Locus | Référence → alternate | Illumina | ONT BQ7 | Fenêtre initiale |
|---|---|---:|---:|---|
| III:105658 | A → T | 0 T / 40 fragments | 7 T / 10 lectures | `III:105558-105758` |
| IV:308249 | T → C | 0 C / 52 fragments | 12 C / 13 lectures | `IV:308149-308349` |

Ces différences de support sont mesurées dans les alignements. Elles ne permettent pas encore de déclarer le SNP ONT vrai ou faux.

Pour chacun, procéder ainsi :

1. Zoomer sur la base exacte et vérifier l'identité de la base dans la référence.
2. Sur chaque piste, examiner les bases alternatives dans les deux orientations, leurs qualités de base et de mapping, les insertions/délétions proches, les extrémités de lectures et le clipping local. Cliquer une lecture pour consulter ses détails.
3. Élargir progressivement à ±500 puis ±2 000 bases pour voir le contexte. La présence d'un clipping quelque part sur une longue lecture n'implique pas un problème au SNP.
4. Noter la version d'IGV, le seuil MAPQ, les réglages des qualités de base, les alignements exclus et le sous-échantillonnage. Pour examiner la distribution des allèles, trier les lectures à la position et colorer par orientation si utile ; relever les réglages choisis. Les deux technologies doivent être visibles sur la capture.
5. Sauvegarder une capture nommée `III_105658.png` ou `IV_308249.png` et écrire une observation descriptive. Exemple : « T visible sur les deux orientations ONT, pas de T dans les fragments Illumina affichés, insertion voisine à examiner ». Ne pas écrire « faux positif confirmé » sans preuve indépendante.

## Renseigner la revue

Copier `evidence/visual_review_template.tsv` sous `visual_review_completed.tsv`. Renseigner `reviewer`, `review_date`, `igv_version`, `notes` et `snapshot_filename`. Dans `visual_status`, utiliser `reviewed` uniquement après examen réel, ou `inconclusive` si l'affichage ne permet pas de conclure. Les autres loci restent `not_reviewed`.

Examiner ensuite **I:2533 G>T** : Illumina 8/14 observations alternatives contre ONT 8/8, avec génotypes archivés 0/1 et 1/1. Puis les trois exemples ONT-only ayant aussi du support Illumina : I:27486, IV:1364942 et IX:37309. Ils illustrent pourquoi l'absence dans un VCF filtré ne signifie pas l'absence d'allèle dans les lectures.

Ajouter la fiche et les captures au dépôt seulement après revue, en conservant les méthodes et les limites. Les observations doivent être attribuées à la personne qui les a effectuées.

Documentation primaire : [sessions IGV](https://igv.org/doc/desktop/UserGuide/sessions/), [SAMtools view et son itérateur multi-région](https://www.htslib.org/doc/samtools-view.html).
