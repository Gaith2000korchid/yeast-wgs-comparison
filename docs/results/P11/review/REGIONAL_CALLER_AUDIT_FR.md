# Audit régional des trois SNPs avec support partagé

Le 1 octobre 2026, une reconstruction ciblée avec **bcftools 1.21 embarqué dans pysam 0.23.3 / HTSlib 1.21** a suivi le pileup, l'appel, la normalisation et le filtrage aux trois sites classés ONT-only malgré leur support Illumina dans la [revue IGV Web](IGV_WEB_REVIEW_FR.md).

**Résultat de cette reconstruction régionale : les trois SNPs Illumina sont appelés et restent inchangés après normalisation, puis sont retirés au filtrage.** Une expérience contrôlée supprimant uniquement BAQ en Illumina fait franchir les filtres à deux de ces trois sites. Il s'agit d'un diagnostic de la méthode, pas d'une recommandation de désactiver BAQ ni d'une validation biologique.

## Où les SNPs disparaissent

Configuration reprise de `config/P11.yaml` : MAPQ≥20, BQ≥13, cap pileup=1000, exclusion des flags 3844, appel diploïde. L'Illumina conserve les réglages BAQ par défaut. La normalisation est `norm -f reference.fa -c e -m -any`. Le filtre retire les records si QUAL manque ou est <30, ou si **FORMAT/DP** manque, est <8 ou >200.

| SNP Illumina | Présent après appel / normalisation | QUAL | INFO/DP brut | FORMAT/DP utilisé par le filtre | AD réf,alt | Motifs de retrait |
|---|---|---:|---:|---:|---|---|
| I:27486 G>A | oui / oui | 59,4147 | 4 | 3 | 0,3 | FORMAT/DP<8 |
| IV:1364942 C>T | oui / oui | 23,434 | 29 | 2 | 0,2 | QUAL<30 et FORMAT/DP<8 |
| IX:37309 A>G | oui / oui | 84,415 | 62 | 6 | 0,6 | FORMAT/DP<8 |

Les trois records portent GT=1/1 dans cette reconstruction. **INFO/DP, FORMAT/DP, profondeur IGV et comptage diagnostique de fragments sont des mesures différentes.** Une profondeur brute élevée ne garantit donc pas de franchir un filtre portant sur les bases retenues par le caller.

Pour contrôle, l'expérience ONT avec BQ≥7, `-B --skip-indels`, conserve les trois SNPs après appel, normalisation et filtrage : FORMAT/DP=15, 8 et 13 ; QUAL=50,208, 42,4147 et 61,4147 ; GT=1/1 aux trois loci. Leur présence ONT et absence Illumina après filtrage reproduisent le classement archivé pour ces trois sites. Les masques de profondeur globaux ne sont pas recalculés ici.

## Expérience contrôlée sur BAQ

Une troisième reconstruction utilise exactement les mêmes BAMs régionaux, référence, fenêtres, seuils et appel diploïde que l'Illumina précédente, en ajoutant **uniquement `-B` à mpileup**. L'appel d'indels reste actif en Illumina dans les deux expériences.

| Locus | FORMAT/DP par défaut → sans BAQ | QUAL par défaut → sans BAQ | Après filtrage sans BAQ |
|---|---:|---:|---|
| I:27486 | 3 → 3 | 59,4147 → 73,4149 | retiré : profondeur <8 |
| IV:1364942 | 2 → 26 | 23,434 → 228,358 | conservé |
| IX:37309 | 6 → 46 | 84,415 → 225,417 | conservé |

BAQ estime l'incertitude de placement des bases dans un alignement et peut réduire leur qualité effective. Le mode par défaut de mpileup utilise des heuristiques pour l'appliquer dans certains contextes. L'ajout de `-B` désactive cet ajustement. Dans cette expérience, il change fortement le nombre de bases retenues à IV:1364942 et IX:37309. Cela démontre une sensibilité au traitement BAQ dans ces fenêtres ; cela ne démontre pas que les bases supplémentaires sont correctes. Les insertions/délétions visibles dans IGV motivent l'examen du contexte, sans prouver la cause locale exacte de l'ajustement.

## Méthode et limites

- Entrées : BAMs régionaux du paquet P11, référence complète, fenêtres d'extraction ±2000 bases conservées à l'identique. Aucun nouvel alignement n'est réalisé.
- Le script appelle les sous-commandes bcftools via pysam, avec des VCF intermédiaires explicites, puis conserve les records aux trois positions, leurs INFO/FORMAT et les commandes. Il vérifie que les motifs de rejet expliquent la présence ou l'absence dans la sortie réelle de `bcftools view`.
- Le [JSON de l'audit](regional_caller_audit.json) conserve les empreintes des BAMs, index, FASTA, FAI, BED et du script, les versions, les paramètres, les messages et les lignes VCF sélectionnées. Les chemins d'exécution sont des informations de provenance, pas des fichiers à télécharger.
- **Ce ne sont pas les intermédiaires originaux du run génome entier.** Le paquet garantit l'identité des comptes diagnostiques aux loci, pas l'identité de tout le comportement du caller. Des mates hors des fenêtres peuvent manquer et le comportement heuristique BAQ peut dépendre du contexte de pileup. Ne pas présenter ces valeurs régionales comme les valeurs historiques exactes du run complet.
- Les valeurs QUAL peuvent dépendre du contexte de régions ; garder le BED fourni pour répéter l'expérience. Une extraction sur la seule base cible n'est pas équivalente à ces fenêtres.
- Nous n'avons pas reconstruit le génome entier, estimé la précision, validé les génotypes ni modifié les paramètres du workflow principal. BAQ désactivé est un contrôle expérimental ; changer les paramètres finaux demanderait une évaluation plus large et une vérité indépendante.

## Reproduire

Extraire le [paquet P11](../../../IGV_REVIEW_FR.md), puis depuis la racine du dépôt :

```bash
micromamba create -y -f environment-review.yml
micromamba run -n yeast-wgs-review python workflow/scripts/audit_regional_calls.py \
  --bundle /chemin/vers/igv_bundle \
  --output results/P11/regional_caller_audit
```

La destination doit être nouvelle. Elle reçoit les VCF intermédiaires locaux et `audit.json`. La version bcftools mesurée figure dans l'en-tête archivé du pileup. Les VCFs de fenêtres complètes ne sont pas ajoutés à Git ; le JSON ne conserve que les positions étudiées.

Un test d'intégration crée deux BAMs indexés synthétiques, appelle deux SNPs et vérifie qu'un SNP de profondeur 3 est présent avant filtrage puis retiré par le seuil FORMAT/DP, tandis qu'un SNP de profondeur 12 est conservé. Les **21 tests unitaires et d'intégration Python** passent localement.

Documentation primaire : [bcftools mpileup, BAQ et filtres](https://samtools.github.io/bcftools/bcftools.html), [pysam : sous-commandes samtools/bcftools](https://pysam.readthedocs.io/en/v0.23.3/usage.html#using-samtools-and-bcftools-commands-within-python).
