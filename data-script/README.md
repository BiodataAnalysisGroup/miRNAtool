Maria's version

To do:

- the script is based on one plate (there are three different). The qc and normalization part will be performed on each plate. The different plates will be merged before differentail expression analysis.

- exclude samples based on: NAs, 'Fail' criterion, user's choise

- the user could select the endogenous or exogenous normalazitaion

- line 212 again with the multimir_results.predicted@summary$target_symbol as input

- the user could select the criterion of GO_enrich and KEGG_enrich


- tables as csv on the working directory:
1. normalized_data
2. data.all (statistical significant miRs)
3. multimir_results.validated@summary$target_symbol and multimir_results.validated@summary$target_symbol as one csv (two different rows)
4. KEGG_enrich.f (for both multimir_results.validated@summary$target_symbol and multimir_results.validated@summary$target_symbol)
5. GO_enrich.f (for both multimir_results.validated@summary$target_symbol and multimir_results.validated@summary$target_symbol)

