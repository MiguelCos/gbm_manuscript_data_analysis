GBM Recurrence - Cox Proportional Hazards Model
================
Justin Sing - Miguel Cosenza-Contreras

- <a href="#initial-data-loading-and-wrangling"
  id="toc-initial-data-loading-and-wrangling"><span
  class="toc-section-number">1</span> Initial data loading and
  wrangling</a>
- <a href="#survival-analysis" id="toc-survival-analysis"><span
  class="toc-section-number">2</span> Survival analysis</a>
  - <a href="#cox-phm---including-only-sex-original-analysis"
    id="toc-cox-phm---including-only-sex-original-analysis"><span
    class="toc-section-number">2.1</span> Cox-PHM - including only sex
    (original analysis)</a>
    - <a href="#genarate-tabular-and-graphical-summaries"
      id="toc-genarate-tabular-and-graphical-summaries"><span
      class="toc-section-number">2.1.1</span> Genarate tabular and graphical
      summaries</a>
  - <a
    href="#cox-phm---including-only-sex-and-age-revision-version-analysis"
    id="toc-cox-phm---including-only-sex-and-age-revision-version-analysis"><span
    class="toc-section-number">2.2</span> Cox-PHM - including only sex and
    age (revision version analysis)</a>

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE)

## Required packages ----
library(tidyverse)
```

    Warning: package 'ggplot2' was built under R version 4.3.1

    Warning: package 'purrr' was built under R version 4.3.1

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)
library(ggrepel)
library(cowplot)
```


    Attaching package: 'cowplot'

    The following object is masked from 'package:ggpubr':

        get_legend

    The following object is masked from 'package:lubridate':

        stamp

``` r
library(here)
```

    here() starts at C:/Users/migue/OneDrive/Documentos/R_Projects/5_projects/gbm/gbm_manuscript_data_analysis

``` r
library(survival)
```

    Warning: package 'survival' was built under R version 4.3.1

``` r
library(gt)


theme_set(theme(axis.text.x = element_text(hjust = 0.5, 
                                           vjust = 0, 
                                           size = 12, 
                                           angle = 360),
                axis.text.y = element_text(hjust = 0.5, 
                                           vjust = 0, 
                                           size = 12),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", 
                                            fill = NA, 
                                            linewidth = 1.5),
                axis.title=element_text(size = 12),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                legend.key.height = unit(4, 
                                         'mm'),
                legend.key.width = unit(4, 
                                        'mm'),
                legend.key.size = unit(10,
                                        'mm'),
                legend.position = "bottom"))
```

# Initial data loading and wrangling

Load the annotation with ASAH1 expression data

``` r
ttr_surv_data <- read_tsv(here::here("data/elisa_plasma_proteomics/ttr_survival_data.tsv"))
```

# Survival analysis

## Cox-PHM - including only sex (original analysis)

Fitting the model:

``` r
## COXPH ----
#ttr_surv_data
fit <- coxph(Surv(Time.to.reccurence..days., status) ~ sex +  ASAH1.ng.ml + SYNM.ng.ml + GPNMB.ng.ml + MMP.9.ng.ml, data = ttr_surv_data)
```

Model summary:

``` r
print(summary(fit))
```

    Call:
    coxph(formula = Surv(Time.to.reccurence..days., status) ~ sex + 
        ASAH1.ng.ml + SYNM.ng.ml + GPNMB.ng.ml + MMP.9.ng.ml, data = ttr_surv_data)

      n= 19, number of events= 18 
       (9 observations deleted due to missingness)

                    coef exp(coef) se(coef)      z Pr(>|z|)  
    sex         -0.03249   0.96803  0.58313 -0.056   0.9556  
    ASAH1.ng.ml -0.52807   0.58974  0.22298 -2.368   0.0179 *
    SYNM.ng.ml   0.06776   1.07011  0.18903  0.358   0.7200  
    GPNMB.ng.ml -0.32392   0.72331  0.33816 -0.958   0.3381  
    MMP.9.ng.ml -0.41722   0.65888  0.27693 -1.507   0.1319  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                exp(coef) exp(-coef) lower .95 upper .95
    sex            0.9680     1.0330    0.3087     3.036
    ASAH1.ng.ml    0.5897     1.6957    0.3809     0.913
    SYNM.ng.ml     1.0701     0.9345    0.7388     1.550
    GPNMB.ng.ml    0.7233     1.3825    0.3728     1.403
    MMP.9.ng.ml    0.6589     1.5177    0.3829     1.134

    Concordance= 0.718  (se = 0.077 )
    Likelihood ratio test= 7.88  on 5 df,   p=0.2
    Wald test            = 7.45  on 5 df,   p=0.2
    Score (logrank) test = 8.17  on 5 df,   p=0.1

Observation: Sex does not seem to represent a significant variable in
the model.

### Genarate tabular and graphical summaries

``` r
library(finalfit)

dependent_var <- "Surv(Time.to.reccurence..days., status)"
explanatory_vars <- c("ASAH1.ng.ml", "SYNM.ng.ml", "GPNMB.ng.ml", "MMP.9.ng.ml")

ttr_surv_data %>% 
    coxphmulti(dependent_var, explanatory_vars) %>% 
    cox.zph() %>% 
    {zph_result <<- .} %>% 
    plot(var=3)
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-5-1.png)

``` r
#zph_result

ttr_surv_data %>%
  finalfit(dependent_var, explanatory_vars) %>%
  gt()
```

<div id="ptzfmqoxzp" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ptzfmqoxzp table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ptzfmqoxzp thead, #ptzfmqoxzp tbody, #ptzfmqoxzp tfoot, #ptzfmqoxzp tr, #ptzfmqoxzp td, #ptzfmqoxzp th {
  border-style: none;
}

#ptzfmqoxzp p {
  margin: 0;
  padding: 0;
}

#ptzfmqoxzp .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ptzfmqoxzp .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ptzfmqoxzp .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ptzfmqoxzp .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ptzfmqoxzp .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ptzfmqoxzp .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ptzfmqoxzp .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ptzfmqoxzp .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ptzfmqoxzp .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ptzfmqoxzp .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ptzfmqoxzp .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ptzfmqoxzp .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ptzfmqoxzp .gt_spanner_row {
  border-bottom-style: hidden;
}

#ptzfmqoxzp .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ptzfmqoxzp .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ptzfmqoxzp .gt_from_md > :first-child {
  margin-top: 0;
}

#ptzfmqoxzp .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ptzfmqoxzp .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ptzfmqoxzp .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ptzfmqoxzp .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ptzfmqoxzp .gt_row_group_first td {
  border-top-width: 2px;
}

#ptzfmqoxzp .gt_row_group_first th {
  border-top-width: 2px;
}

#ptzfmqoxzp .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ptzfmqoxzp .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ptzfmqoxzp .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ptzfmqoxzp .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ptzfmqoxzp .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ptzfmqoxzp .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ptzfmqoxzp .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ptzfmqoxzp .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ptzfmqoxzp .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ptzfmqoxzp .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ptzfmqoxzp .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ptzfmqoxzp .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ptzfmqoxzp .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ptzfmqoxzp .gt_left {
  text-align: left;
}

#ptzfmqoxzp .gt_center {
  text-align: center;
}

#ptzfmqoxzp .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ptzfmqoxzp .gt_font_normal {
  font-weight: normal;
}

#ptzfmqoxzp .gt_font_bold {
  font-weight: bold;
}

#ptzfmqoxzp .gt_font_italic {
  font-style: italic;
}

#ptzfmqoxzp .gt_super {
  font-size: 65%;
}

#ptzfmqoxzp .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ptzfmqoxzp .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ptzfmqoxzp .gt_indent_1 {
  text-indent: 5px;
}

#ptzfmqoxzp .gt_indent_2 {
  text-indent: 10px;
}

#ptzfmqoxzp .gt_indent_3 {
  text-indent: 15px;
}

#ptzfmqoxzp .gt_indent_4 {
  text-indent: 20px;
}

#ptzfmqoxzp .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Dependent: Surv(Time.to.reccurence..days., status)">Dependent: Surv(Time.to.reccurence..days., status)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=" "> </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="all">all</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="HR (univariable)">HR (univariable)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="HR (multivariable)">HR (multivariable)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">ASAH1.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">2.2 (1.5)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.77 (0.58-1.02, p=0.066)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.59 (0.39-0.91, p=0.016)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">SYNM.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">9.8 (2.2)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.91 (0.79-1.05, p=0.183)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">1.07 (0.75-1.54, p=0.703)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">GPNMB.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">5.3 (1.5)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.85 (0.65-1.10, p=0.213)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.72 (0.37-1.39, p=0.331)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">MMP.9.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">7.5 (1.2)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.87 (0.57-1.34, p=0.536)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.65 (0.40-1.08, p=0.097)</td></tr>
  </tbody>
  
  
</table>
</div>

``` r
ttr_surv_data %>%
  hr_plot(dependent_var, explanatory_vars, 
          dependent_label="Time to Recurrence")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-5-2.png)

``` r
ttr_surv_data %>%
  ff_plot(dependent_var, explanatory_vars, 
          dependent_label="Time to Recurrence")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-5-3.png)

## Cox-PHM - including only sex and age (revision version analysis)

This section was created as part of the revision process, to check for
the effect of age and sex in the model.

``` r
## COXPH ----
#ttr_surv_data
fit2 <- coxph(Surv(Time.to.reccurence..days., status) ~ sex + Age.at.surgery..years. + ASAH1.ng.ml + SYNM.ng.ml + GPNMB.ng.ml + MMP.9.ng.ml, data = ttr_surv_data)
```

Model summary:

``` r
print(summary(fit2))
```

    Call:
    coxph(formula = Surv(Time.to.reccurence..days., status) ~ sex + 
        Age.at.surgery..years. + ASAH1.ng.ml + SYNM.ng.ml + GPNMB.ng.ml + 
        MMP.9.ng.ml, data = ttr_surv_data)

      n= 19, number of events= 18 
       (9 observations deleted due to missingness)

                               coef exp(coef) se(coef)      z Pr(>|z|)  
    sex                    -0.22790   0.79620  0.62276 -0.366   0.7144  
    Age.at.surgery..years.  0.06170   1.06365  0.03614  1.707   0.0878 .
    ASAH1.ng.ml            -0.60744   0.54474  0.24005 -2.530   0.0114 *
    SYNM.ng.ml              0.12710   1.13553  0.18555  0.685   0.4934  
    GPNMB.ng.ml            -0.26025   0.77086  0.32028 -0.813   0.4165  
    MMP.9.ng.ml            -0.39958   0.67060  0.29407 -1.359   0.1742  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                           exp(coef) exp(-coef) lower .95 upper .95
    sex                       0.7962     1.2560    0.2349     2.698
    Age.at.surgery..years.    1.0636     0.9402    0.9909     1.142
    ASAH1.ng.ml               0.5447     1.8357    0.3403     0.872
    SYNM.ng.ml                1.1355     0.8806    0.7893     1.634
    GPNMB.ng.ml               0.7709     1.2973    0.4115     1.444
    MMP.9.ng.ml               0.6706     1.4912    0.3768     1.193

    Concordance= 0.735  (se = 0.076 )
    Likelihood ratio test= 10.88  on 6 df,   p=0.09
    Wald test            = 8.77  on 6 df,   p=0.2
    Score (logrank) test = 10.3  on 6 df,   p=0.1

``` r
?coxphmulti
```

Observation: p-value for age is 0.0878 in the model; while ASAH1
concentration is stil below 0.05.

``` r
dependent_var <- "Surv(Time.to.reccurence..days., status)"
explanatory_vars2 <- c("ASAH1.ng.ml", "SYNM.ng.ml", "GPNMB.ng.ml", "MMP.9.ng.ml",
                      "sex", "Age.at.surgery..years.") # includes sex and age

ttr_surv_data %>% 
    coxphmulti(dependent_var, explanatory_vars2) %>% 
    cox.zph() %>% 
    {zph_result2 <<- .} %>% 
    plot(var=3)
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-9-1.png)

``` r
#zph_result2

ttr_surv_data %>%
  finalfit(dependent_var, explanatory_vars2) %>%
  gt()
```

<div id="usyelbpbzf" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#usyelbpbzf table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#usyelbpbzf thead, #usyelbpbzf tbody, #usyelbpbzf tfoot, #usyelbpbzf tr, #usyelbpbzf td, #usyelbpbzf th {
  border-style: none;
}

#usyelbpbzf p {
  margin: 0;
  padding: 0;
}

#usyelbpbzf .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#usyelbpbzf .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#usyelbpbzf .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#usyelbpbzf .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#usyelbpbzf .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#usyelbpbzf .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#usyelbpbzf .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#usyelbpbzf .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#usyelbpbzf .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#usyelbpbzf .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#usyelbpbzf .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#usyelbpbzf .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#usyelbpbzf .gt_spanner_row {
  border-bottom-style: hidden;
}

#usyelbpbzf .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#usyelbpbzf .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#usyelbpbzf .gt_from_md > :first-child {
  margin-top: 0;
}

#usyelbpbzf .gt_from_md > :last-child {
  margin-bottom: 0;
}

#usyelbpbzf .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#usyelbpbzf .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#usyelbpbzf .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#usyelbpbzf .gt_row_group_first td {
  border-top-width: 2px;
}

#usyelbpbzf .gt_row_group_first th {
  border-top-width: 2px;
}

#usyelbpbzf .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#usyelbpbzf .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#usyelbpbzf .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#usyelbpbzf .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#usyelbpbzf .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#usyelbpbzf .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#usyelbpbzf .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#usyelbpbzf .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#usyelbpbzf .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#usyelbpbzf .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#usyelbpbzf .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#usyelbpbzf .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#usyelbpbzf .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#usyelbpbzf .gt_left {
  text-align: left;
}

#usyelbpbzf .gt_center {
  text-align: center;
}

#usyelbpbzf .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#usyelbpbzf .gt_font_normal {
  font-weight: normal;
}

#usyelbpbzf .gt_font_bold {
  font-weight: bold;
}

#usyelbpbzf .gt_font_italic {
  font-style: italic;
}

#usyelbpbzf .gt_super {
  font-size: 65%;
}

#usyelbpbzf .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#usyelbpbzf .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#usyelbpbzf .gt_indent_1 {
  text-indent: 5px;
}

#usyelbpbzf .gt_indent_2 {
  text-indent: 10px;
}

#usyelbpbzf .gt_indent_3 {
  text-indent: 15px;
}

#usyelbpbzf .gt_indent_4 {
  text-indent: 20px;
}

#usyelbpbzf .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Dependent: Surv(Time.to.reccurence..days., status)">Dependent: Surv(Time.to.reccurence..days., status)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=" "> </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="all">all</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="HR (univariable)">HR (univariable)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="HR (multivariable)">HR (multivariable)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">ASAH1.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">2.2 (1.5)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.77 (0.58-1.02, p=0.066)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.54 (0.34-0.87, p=0.011)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">SYNM.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">9.8 (2.2)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.91 (0.79-1.05, p=0.183)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">1.14 (0.79-1.63, p=0.493)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">GPNMB.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">5.3 (1.5)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.85 (0.65-1.10, p=0.213)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.77 (0.41-1.44, p=0.416)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">MMP.9.ng.ml</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">7.5 (1.2)</td>
<td headers="HR (univariable)" class="gt_row gt_left">0.87 (0.57-1.34, p=0.536)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.67 (0.38-1.19, p=0.174)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">sex</td>
<td headers="" class="gt_row gt_left">1</td>
<td headers="all" class="gt_row gt_right">12 (42.9)</td>
<td headers="HR (univariable)" class="gt_row gt_left">-</td>
<td headers="HR (multivariable)" class="gt_row gt_left">-</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left"></td>
<td headers="" class="gt_row gt_left">2</td>
<td headers="all" class="gt_row gt_right">16 (57.1)</td>
<td headers="HR (univariable)" class="gt_row gt_left">1.35 (0.56-3.22, p=0.504)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">0.80 (0.23-2.70, p=0.714)</td></tr>
    <tr><td headers="Dependent: Surv(Time.to.reccurence..days., status)" class="gt_row gt_left">Age.at.surgery..years.</td>
<td headers="" class="gt_row gt_left">Mean (SD)</td>
<td headers="all" class="gt_row gt_right">62.4 (8.8)</td>
<td headers="HR (univariable)" class="gt_row gt_left">1.06 (1.01-1.12, p=0.026)</td>
<td headers="HR (multivariable)" class="gt_row gt_left">1.06 (0.99-1.14, p=0.088)</td></tr>
  </tbody>
  
  
</table>
</div>

``` r
ttr_surv_data %>%
  hr_plot(dependent_var, explanatory_vars2, dependent_label="Time to Recurrence")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-9-2.png)

``` r
ttr_surv_data %>%
  ff_plot(dependent_var, explanatory_vars2, dependent_label="Time to Recurrence")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-9-3.png)

When including Age and Sex in the model, ASAH1 concentrations are still
significantly associated with the time to recurrence. Age and sex are
not significant in the model.

``` r
ggsave(here::here("figures/cox_phm_summary_rev_version_w_age_n_sex.pdf"), 
       ttr_surv_data %>%
  ff_plot(dependent_var, explanatory_vars2, dependent_label="Time to Recurrence"), 
       width = 20, 
       height = 20,
       units = "cm")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-10-1.png)

``` r
ggsave(here::here("figures/cox_phm_summary_rev_version_w_age_n_sex.png"), 
       ttr_surv_data %>%
  ff_plot(dependent_var, explanatory_vars2, dependent_label="Time to Recurrence"), 
       width = 20, 
       height = 20,
       units = "cm")
```

![](gbm_recurrence_cox_phm_incl_age_sex_revision_files/figure-gfm/unnamed-chunk-11-1.png)

``` r
dev.off()
```

    null device 
              1 
