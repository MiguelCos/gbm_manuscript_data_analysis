GBM Proteomics - Patient summary
================
Miguel Cosenza-Contreras
2023-02-16

- <a href="#summary-table-of-patients-minimal-information"
  id="toc-summary-table-of-patients-minimal-information"><span
  class="toc-section-number">1</span> Summary table of patients (minimal
  information)</a>
- <a href="#summary-table-of-patients-all-clinical-information"
  id="toc-summary-table-of-patients-all-clinical-information"><span
  class="toc-section-number">2</span> Summary table of patients (all
  clinical information)</a>
- <a href="#summary-table-of-patients-for-proteomics"
  id="toc-summary-table-of-patients-for-proteomics"><span
  class="toc-section-number">3</span> Summary table of patients for
  proteomics</a>

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE)
```

# Required packages

``` r
library(readxl)
library(gtsummary)
library(here)
library(tidyverse)
library(flextable)
```

# Load and prep the data

``` r
patients_proteomics <- patient_table %>%
  filter(`Mass Spectrometry` == "X") %>%
  dplyr::select(Sex, 
                `Age at surgery (years)`, 
                `Time-to-reccurence (days)`,
                `Type of resection`,
                `Tumor localization`, 
                Treatment,
                sample_id = `LAB-ID`) %>%
  mutate(`Time-to-reccurence (days)` = as.numeric(`Time-to-reccurence (days)`))
```

``` r
patients_clinical <- patient_table %>%
  mutate(`Time-to-reccurence (days)` = as.numeric(`Time-to-reccurence (days)`),
         `Age at surgery (years)` = as.numeric(`Age at surgery (years)`),
         `Latency (days)` = as.numeric(`Latency (days)`)) %>% 
  pivot_longer(cols = c("Mass Spectrometry",
                        "Lipidomics",
                        "RT-qPCR",
                        "IHC",
                        "ELISA"),
               names_to = "Experiment",
               values_to = "X") %>%
  filter(!is.na(X)) %>%
  dplyr::select(-c(`LAB-ID`, Number))
```

``` r
patients_clinical_min <- patients_clinical %>%
  dplyr::select(Sex, 
                `Age at surgery (years)`, 
                `Time-to-reccurence (days)`,
                `Latency (days)`,
                `Type of resection`,
                Experiment)
```

# Summary table of patients (minimal information)

``` r
clin_tb_min <- tbl_summary(data = patients_clinical_min,
                           by = Experiment) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 
```

``` r
print(clin_tb_min)
```

    <div id="qabfrwfciq" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>html {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    }

    #qabfrwfciq .gt_table {
      display: table;
      border-collapse: collapse;
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

    #qabfrwfciq .gt_heading {
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

    #qabfrwfciq .gt_title {
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

    #qabfrwfciq .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 0;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #qabfrwfciq .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #qabfrwfciq .gt_col_headings {
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

    #qabfrwfciq .gt_col_heading {
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

    #qabfrwfciq .gt_column_spanner_outer {
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

    #qabfrwfciq .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #qabfrwfciq .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #qabfrwfciq .gt_column_spanner {
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

    #qabfrwfciq .gt_group_heading {
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
    }

    #qabfrwfciq .gt_empty_group_heading {
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

    #qabfrwfciq .gt_from_md > :first-child {
      margin-top: 0;
    }

    #qabfrwfciq .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #qabfrwfciq .gt_row {
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

    #qabfrwfciq .gt_stub {
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

    #qabfrwfciq .gt_stub_row_group {
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

    #qabfrwfciq .gt_row_group_first td {
      border-top-width: 2px;
    }

    #qabfrwfciq .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #qabfrwfciq .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #qabfrwfciq .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #qabfrwfciq .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #qabfrwfciq .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #qabfrwfciq .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #qabfrwfciq .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #qabfrwfciq .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #qabfrwfciq .gt_footnotes {
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

    #qabfrwfciq .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-left: 4px;
      padding-right: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #qabfrwfciq .gt_sourcenotes {
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

    #qabfrwfciq .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #qabfrwfciq .gt_left {
      text-align: left;
    }

    #qabfrwfciq .gt_center {
      text-align: center;
    }

    #qabfrwfciq .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #qabfrwfciq .gt_font_normal {
      font-weight: normal;
    }

    #qabfrwfciq .gt_font_bold {
      font-weight: bold;
    }

    #qabfrwfciq .gt_font_italic {
      font-style: italic;
    }

    #qabfrwfciq .gt_super {
      font-size: 65%;
    }

    #qabfrwfciq .gt_footnote_marks {
      font-style: italic;
      font-weight: normal;
      font-size: 75%;
      vertical-align: 0.4em;
    }

    #qabfrwfciq .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #qabfrwfciq .gt_indent_1 {
      text-indent: 5px;
    }

    #qabfrwfciq .gt_indent_2 {
      text-indent: 10px;
    }

    #qabfrwfciq .gt_indent_3 {
      text-indent: 15px;
    }

    #qabfrwfciq .gt_indent_4 {
      text-indent: 20px;
    }

    #qabfrwfciq .gt_indent_5 {
      text-indent: 25px;
    }
    </style>
      <table class="gt_table">
      
      <thead class="gt_col_headings">
        <tr>
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Variable</strong></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>N</strong></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>ELISA</strong>, N = 30<sup class="gt_footnote_marks">1</sup></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>IHC</strong>, N = 10<sup class="gt_footnote_marks">1</sup></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Lipidomics</strong>, N = 10<sup class="gt_footnote_marks">1</sup></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Mass Spectrometry</strong>, N = 11<sup class="gt_footnote_marks">1</sup></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>RT-qPCR</strong>, N = 20<sup class="gt_footnote_marks">1</sup></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td class="gt_row gt_left" style="font-weight: bold;">Sex</td>
    <td class="gt_row gt_center">81</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td></tr>
        <tr><td class="gt_row gt_left">    m</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">14 (47%)</td>
    <td class="gt_row gt_center">4 (40%)</td>
    <td class="gt_row gt_center">5 (50%)</td>
    <td class="gt_row gt_center">6 (55%)</td>
    <td class="gt_row gt_center">11 (55%)</td></tr>
        <tr><td class="gt_row gt_left">    w</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">16 (53%)</td>
    <td class="gt_row gt_center">6 (60%)</td>
    <td class="gt_row gt_center">5 (50%)</td>
    <td class="gt_row gt_center">5 (45%)</td>
    <td class="gt_row gt_center">9 (45%)</td></tr>
        <tr><td class="gt_row gt_left" style="font-weight: bold;">Age at surgery (years)</td>
    <td class="gt_row gt_center">81</td>
    <td class="gt_row gt_center">62 (54, 68)</td>
    <td class="gt_row gt_center">51 (50, 61)</td>
    <td class="gt_row gt_center">62 (54, 63)</td>
    <td class="gt_row gt_center">51 (46, 58)</td>
    <td class="gt_row gt_center">61 (52, 64)</td></tr>
        <tr><td class="gt_row gt_left" style="font-weight: bold;">Time-to-reccurence (days)</td>
    <td class="gt_row gt_center">74</td>
    <td class="gt_row gt_center">224 (166, 288)</td>
    <td class="gt_row gt_center">284 (160, 456)</td>
    <td class="gt_row gt_center">218 (190, 273)</td>
    <td class="gt_row gt_center">237 (182, 398)</td>
    <td class="gt_row gt_center">288 (192, 418)</td></tr>
        <tr><td class="gt_row gt_left">    Unknown</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">7</td>
    <td class="gt_row gt_center">0</td>
    <td class="gt_row gt_center">0</td>
    <td class="gt_row gt_center">0</td>
    <td class="gt_row gt_center">0</td></tr>
        <tr><td class="gt_row gt_left" style="font-weight: bold;">Latency (days)</td>
    <td class="gt_row gt_center">47</td>
    <td class="gt_row gt_center">437 (336, 603)</td>
    <td class="gt_row gt_center">576 (410, 775)</td>
    <td class="gt_row gt_center">546 (431, 766)</td>
    <td class="gt_row gt_center">681 (516, 1,073)</td>
    <td class="gt_row gt_center">602 (491, 1,054)</td></tr>
        <tr><td class="gt_row gt_left">    Unknown</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">11</td>
    <td class="gt_row gt_center">3</td>
    <td class="gt_row gt_center">6</td>
    <td class="gt_row gt_center">4</td>
    <td class="gt_row gt_center">10</td></tr>
        <tr><td class="gt_row gt_left" style="font-weight: bold;">Type of resection</td>
    <td class="gt_row gt_center">81</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center"></td></tr>
        <tr><td class="gt_row gt_left">    subtotal</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">12 (40%)</td>
    <td class="gt_row gt_center">4 (40%)</td>
    <td class="gt_row gt_center">2 (20%)</td>
    <td class="gt_row gt_center">5 (45%)</td>
    <td class="gt_row gt_center">9 (45%)</td></tr>
        <tr><td class="gt_row gt_left">    total</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">18 (60%)</td>
    <td class="gt_row gt_center">5 (50%)</td>
    <td class="gt_row gt_center">8 (80%)</td>
    <td class="gt_row gt_center">6 (55%)</td>
    <td class="gt_row gt_center">11 (55%)</td></tr>
        <tr><td class="gt_row gt_left">    unknown</td>
    <td class="gt_row gt_center"></td>
    <td class="gt_row gt_center">0 (0%)</td>
    <td class="gt_row gt_center">1 (10%)</td>
    <td class="gt_row gt_center">0 (0%)</td>
    <td class="gt_row gt_center">0 (0%)</td>
    <td class="gt_row gt_center">0 (0%)</td></tr>
      </tbody>
      
      <tfoot class="gt_footnotes">
        <tr>
          <td class="gt_footnote" colspan="7"><sup class="gt_footnote_marks">1</sup> n (%); Median (IQR)</td>
        </tr>
      </tfoot>
    </table>
    </div>

# Summary table of patients (all clinical information)

``` r
tbl_summary(data = patients_clinical,
            by = Experiment) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 
```

<div id="uxfzopyxlg" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#uxfzopyxlg .gt_table {
  display: table;
  border-collapse: collapse;
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

#uxfzopyxlg .gt_heading {
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

#uxfzopyxlg .gt_title {
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

#uxfzopyxlg .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#uxfzopyxlg .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uxfzopyxlg .gt_col_headings {
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

#uxfzopyxlg .gt_col_heading {
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

#uxfzopyxlg .gt_column_spanner_outer {
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

#uxfzopyxlg .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#uxfzopyxlg .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#uxfzopyxlg .gt_column_spanner {
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

#uxfzopyxlg .gt_group_heading {
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
}

#uxfzopyxlg .gt_empty_group_heading {
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

#uxfzopyxlg .gt_from_md > :first-child {
  margin-top: 0;
}

#uxfzopyxlg .gt_from_md > :last-child {
  margin-bottom: 0;
}

#uxfzopyxlg .gt_row {
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

#uxfzopyxlg .gt_stub {
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

#uxfzopyxlg .gt_stub_row_group {
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

#uxfzopyxlg .gt_row_group_first td {
  border-top-width: 2px;
}

#uxfzopyxlg .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uxfzopyxlg .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#uxfzopyxlg .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#uxfzopyxlg .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uxfzopyxlg .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uxfzopyxlg .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#uxfzopyxlg .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#uxfzopyxlg .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uxfzopyxlg .gt_footnotes {
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

#uxfzopyxlg .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#uxfzopyxlg .gt_sourcenotes {
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

#uxfzopyxlg .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#uxfzopyxlg .gt_left {
  text-align: left;
}

#uxfzopyxlg .gt_center {
  text-align: center;
}

#uxfzopyxlg .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#uxfzopyxlg .gt_font_normal {
  font-weight: normal;
}

#uxfzopyxlg .gt_font_bold {
  font-weight: bold;
}

#uxfzopyxlg .gt_font_italic {
  font-style: italic;
}

#uxfzopyxlg .gt_super {
  font-size: 65%;
}

#uxfzopyxlg .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#uxfzopyxlg .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#uxfzopyxlg .gt_indent_1 {
  text-indent: 5px;
}

#uxfzopyxlg .gt_indent_2 {
  text-indent: 10px;
}

#uxfzopyxlg .gt_indent_3 {
  text-indent: 15px;
}

#uxfzopyxlg .gt_indent_4 {
  text-indent: 20px;
}

#uxfzopyxlg .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Variable</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>N</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>ELISA</strong>, N = 30<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>IHC</strong>, N = 10<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Lipidomics</strong>, N = 10<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Mass Spectrometry</strong>, N = 11<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>RT-qPCR</strong>, N = 20<sup class="gt_footnote_marks">1</sup></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Age at surgery (years)</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center">62 (54, 68)</td>
<td class="gt_row gt_center">51 (50, 61)</td>
<td class="gt_row gt_center">62 (54, 63)</td>
<td class="gt_row gt_center">51 (46, 58)</td>
<td class="gt_row gt_center">61 (52, 64)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Sex</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    m</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">14 (47%)</td>
<td class="gt_row gt_center">4 (40%)</td>
<td class="gt_row gt_center">5 (50%)</td>
<td class="gt_row gt_center">6 (55%)</td>
<td class="gt_row gt_center">11 (55%)</td></tr>
    <tr><td class="gt_row gt_left">    w</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">16 (53%)</td>
<td class="gt_row gt_center">6 (60%)</td>
<td class="gt_row gt_center">5 (50%)</td>
<td class="gt_row gt_center">5 (45%)</td>
<td class="gt_row gt_center">9 (45%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Tumor localization</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    left-parietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    left frontal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">4 (13%)</td>
<td class="gt_row gt_center">3 (30%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">2 (18%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    left gyrus cinguli</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    left occipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    left parietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    left parieto-occipitale</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    left tempoparietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    left tempopomesial</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    left temporal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">5 (17%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">2 (18%)</td>
<td class="gt_row gt_center">3 (15%)</td></tr>
    <tr><td class="gt_row gt_left">    left temporal insulär</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    left trigonal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    right frontal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">4 (20%)</td></tr>
    <tr><td class="gt_row gt_left">    right okzipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    right parietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    right parieto-occipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    right parietooccipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    right tempoparietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    right temporal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">7 (23%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    right temporo parietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    right temporooccipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    right trigonal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    right zerebral</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Type of resection</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    subtotal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">12 (40%)</td>
<td class="gt_row gt_center">4 (40%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">5 (45%)</td>
<td class="gt_row gt_center">9 (45%)</td></tr>
    <tr><td class="gt_row gt_left">    total</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">18 (60%)</td>
<td class="gt_row gt_center">5 (50%)</td>
<td class="gt_row gt_center">8 (80%)</td>
<td class="gt_row gt_center">6 (55%)</td>
<td class="gt_row gt_center">11 (55%)</td></tr>
    <tr><td class="gt_row gt_left">    unknown</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">MGMT promotor methylation</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    -</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">10 (33%)</td>
<td class="gt_row gt_center">7 (70%)</td>
<td class="gt_row gt_center">6 (60%)</td>
<td class="gt_row gt_center">8 (73%)</td>
<td class="gt_row gt_center">9 (45%)</td></tr>
    <tr><td class="gt_row gt_left">    (+)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">10 (33%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (15%)</td></tr>
    <tr><td class="gt_row gt_left">    +</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">10 (33%)</td>
<td class="gt_row gt_center">3 (30%)</td>
<td class="gt_row gt_center">3 (30%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">8 (40%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">p53 accumulation</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    -</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    (+)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">11 (37%)</td>
<td class="gt_row gt_center">5 (50%)</td>
<td class="gt_row gt_center">7 (70%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">10 (50%)</td></tr>
    <tr><td class="gt_row gt_left">    +</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">14 (47%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">4 (36%)</td>
<td class="gt_row gt_center">7 (35%)</td></tr>
    <tr><td class="gt_row gt_left">    ++</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (10%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    +++</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Ki67-Li</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    0.05</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    0.1</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    0.15</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">5 (17%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (30%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (15%)</td></tr>
    <tr><td class="gt_row gt_left">    0.2</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">8 (27%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">2 (18%)</td>
<td class="gt_row gt_center">5 (25%)</td></tr>
    <tr><td class="gt_row gt_left">    0.25</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    0.3</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">5 (17%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">3 (30%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">6 (30%)</td></tr>
    <tr><td class="gt_row gt_left">    0.4</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">2 (18%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    0.5</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (10%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    0.6</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    0.75</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    0.9</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    15-20%</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    15-30%</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    30-40%</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">EGFR vIII</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    -</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">25 (83%)</td>
<td class="gt_row gt_center">8 (80%)</td>
<td class="gt_row gt_center">7 (70%)</td>
<td class="gt_row gt_center">10 (91%)</td>
<td class="gt_row gt_center">15 (75%)</td></tr>
    <tr><td class="gt_row gt_left">    (+)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    +</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">4 (13%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (15%)</td></tr>
    <tr><td class="gt_row gt_left">    ++</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Treatment</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    /</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    Chemo external</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    keine RT / keine Chemo - bei uns "nur" Rezidiv OP</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    no RT/Chemo only for iGBM</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (27,2 Gy) - Abbruch</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (34 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (36 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (40,5 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (40,5 Gy) /Chemo (12 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (40,5 Gy) /Chemo (3 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (40,5 Gy) /Chemo (6 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (50,5 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">3 (27%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (1 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (12 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (2 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (6.7%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (3 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (10%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">2 (20%)</td>
<td class="gt_row gt_center">1 (9.1%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (4 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">4 (20%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (5 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (6 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">6 (20%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">2 (18%)</td>
<td class="gt_row gt_center">2 (10%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (7 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo ?(1 Zyklus)?</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo nur zur RT</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (3.3%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left">    RT extern</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (10%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (5.0%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Time-to-reccurence (days)</td>
<td class="gt_row gt_center">74</td>
<td class="gt_row gt_center">224 (166, 288)</td>
<td class="gt_row gt_center">284 (160, 456)</td>
<td class="gt_row gt_center">218 (190, 273)</td>
<td class="gt_row gt_center">237 (182, 398)</td>
<td class="gt_row gt_center">288 (192, 418)</td></tr>
    <tr><td class="gt_row gt_left">    Unknown</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">7</td>
<td class="gt_row gt_center">0</td>
<td class="gt_row gt_center">0</td>
<td class="gt_row gt_center">0</td>
<td class="gt_row gt_center">0</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Latency (days)</td>
<td class="gt_row gt_center">47</td>
<td class="gt_row gt_center">437 (336, 603)</td>
<td class="gt_row gt_center">576 (410, 775)</td>
<td class="gt_row gt_center">546 (431, 766)</td>
<td class="gt_row gt_center">681 (516, 1,073)</td>
<td class="gt_row gt_center">602 (491, 1,054)</td></tr>
    <tr><td class="gt_row gt_left">    Unknown</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center">3</td>
<td class="gt_row gt_center">6</td>
<td class="gt_row gt_center">4</td>
<td class="gt_row gt_center">10</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">X</td>
<td class="gt_row gt_center">81</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    X</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">30 (100%)</td>
<td class="gt_row gt_center">10 (100%)</td>
<td class="gt_row gt_center">10 (100%)</td>
<td class="gt_row gt_center">11 (100%)</td>
<td class="gt_row gt_center">20 (100%)</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="7"><sup class="gt_footnote_marks">1</sup> Median (IQR); n (%)</td>
    </tr>
  </tfoot>
</table>
</div>

# Summary table of patients for proteomics

``` r
tbl_summary(data = patients_proteomics) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 
```

<div id="hmqvnieuqh" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#hmqvnieuqh .gt_table {
  display: table;
  border-collapse: collapse;
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

#hmqvnieuqh .gt_heading {
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

#hmqvnieuqh .gt_title {
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

#hmqvnieuqh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#hmqvnieuqh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hmqvnieuqh .gt_col_headings {
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

#hmqvnieuqh .gt_col_heading {
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

#hmqvnieuqh .gt_column_spanner_outer {
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

#hmqvnieuqh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#hmqvnieuqh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#hmqvnieuqh .gt_column_spanner {
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

#hmqvnieuqh .gt_group_heading {
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
}

#hmqvnieuqh .gt_empty_group_heading {
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

#hmqvnieuqh .gt_from_md > :first-child {
  margin-top: 0;
}

#hmqvnieuqh .gt_from_md > :last-child {
  margin-bottom: 0;
}

#hmqvnieuqh .gt_row {
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

#hmqvnieuqh .gt_stub {
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

#hmqvnieuqh .gt_stub_row_group {
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

#hmqvnieuqh .gt_row_group_first td {
  border-top-width: 2px;
}

#hmqvnieuqh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hmqvnieuqh .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#hmqvnieuqh .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#hmqvnieuqh .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hmqvnieuqh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hmqvnieuqh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#hmqvnieuqh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#hmqvnieuqh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hmqvnieuqh .gt_footnotes {
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

#hmqvnieuqh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#hmqvnieuqh .gt_sourcenotes {
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

#hmqvnieuqh .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#hmqvnieuqh .gt_left {
  text-align: left;
}

#hmqvnieuqh .gt_center {
  text-align: center;
}

#hmqvnieuqh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#hmqvnieuqh .gt_font_normal {
  font-weight: normal;
}

#hmqvnieuqh .gt_font_bold {
  font-weight: bold;
}

#hmqvnieuqh .gt_font_italic {
  font-style: italic;
}

#hmqvnieuqh .gt_super {
  font-size: 65%;
}

#hmqvnieuqh .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#hmqvnieuqh .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#hmqvnieuqh .gt_indent_1 {
  text-indent: 5px;
}

#hmqvnieuqh .gt_indent_2 {
  text-indent: 10px;
}

#hmqvnieuqh .gt_indent_3 {
  text-indent: 15px;
}

#hmqvnieuqh .gt_indent_4 {
  text-indent: 20px;
}

#hmqvnieuqh .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Variable</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>N</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>N = 11</strong><sup class="gt_footnote_marks">1</sup></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Sex</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    m</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">6 (55%)</td></tr>
    <tr><td class="gt_row gt_left">    w</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">5 (45%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Age at surgery (years)</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center">51 (46, 58)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Time-to-reccurence (days)</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center">237 (182, 398)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Type of resection</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    subtotal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">5 (45%)</td></tr>
    <tr><td class="gt_row gt_left">    total</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">6 (55%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Tumor localization</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    left frontal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (18%)</td></tr>
    <tr><td class="gt_row gt_left">    left parieto-occipitale</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    left temporal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (18%)</td></tr>
    <tr><td class="gt_row gt_left">    left temporal insulär</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    right okzipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    right parietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    right parieto-occipital</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    right tempoparietal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    right trigonal</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">Treatment</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    Chemo external</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (36 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (27%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">3 (27%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (3 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    RT (60 Gy) /Chemo (6 cycles)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (18%)</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">sample_id</td>
<td class="gt_row gt_center">11</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left">    2012/49</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2013/210</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2013/84</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2015/139</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2015/150</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2015/161</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2015/18</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2016/145</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">2 (18%)</td></tr>
    <tr><td class="gt_row gt_left">    2016/65</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
    <tr><td class="gt_row gt_left">    2016/70</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center">1 (9.1%)</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="3"><sup class="gt_footnote_marks">1</sup> n (%); Median (IQR)</td>
    </tr>
  </tfoot>
</table>
</div>
