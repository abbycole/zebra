---
output: 
  pdf_document: 
    latex_engine: xelatex
mainfont: Calibri
fig_caption: true
params: 
    set_title: "My Title!"
date: "31 July, 2017"
---

---
title: "Microbiome Diet Study Participant Report MCTs01"
---

Thank you for participating in the Knights Lab's citizen science project! 

We have made you a personalized report from your samples and dietary data. Included in this report is information about your daily dietary intake and your daily microbiome variation. This report is not intended to medical advice or to be used to diagnose or treat disease. If you have any questions about your report, or you would like access to your raw sequence data, please contact the study coordinator Abby at cole0463@umn.edu.

#Your Dietary Intake Data

This first table and Figure 1 below show your average macronutrient intake as a percentage of calories and the typical macronutrient intake ranges for US adults:


| Type                | Average Macronutrient Intake                               | Typical Intake Range        |
|---------------------|------------------------------------------------------------|-----------------------------|
| CALORIES (kcal)     | 2233.97                            |   1800 - 3200               |
| PROTEIN (%)         | 25.21%        |   10 - 35%                  | 
| TOTAL FAT (%)       | 17.06%        |   20 - 35%                  |  
| CARBS (%)           | 57.73%        |   45 - 65%                  |                                   

Based on your height, weight, gender, age, and self-reported activity level you need approximately 2100 kcal per day to maintain your current weight. During the study your weight change was -0.4 kg. 

##Figure 1:

\begin{center}\includegraphics{/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-8-1} \end{center}

**Figure 1** shows the day to day variation in your consumption of macronutrients viewed as a percentage of total kilocalorie intake.

Please note, it is common for people who follow some diets, like the low carbohydrate diet, to have intake ranges that differ from the typical intake ranges. The next table will show if your diet is nutritionally adequate.


\newpage
This second table shows your average daily intake of key micronutrients from food. Also shown are the average intake levels for all of our study participants and the recommended dietary allowance levels by gender. If you take a vitamin and mineral supplement, your intake from that supplement is not reflected here:


| Nutrient            | Your Average                    | Study Average                                       |Recommended Dietary Allowance(Male/Female)|
|---------------------|---------------------------------|-----------------------------------------------------|------------------------------------------|          
| FOLATE (mcg)        | 677.35 | 435.99 |                 400                      |
| CHOLINE (mg)        | 415.25| 364.25|                550 / 425                 |
| CALCIUM (mg)        | 1477.48 | 1054.42 |                 1000 *                   |
| SODIUM (mg)         | 4413.8 | 3429.46 |                 1500 *                   |
| POTASSIUM (mg)      | 3768.93 | 2904.59 |                 4700                     | 
| MAGNESIUM (mg)      | 402.16 | 364.25 |               400 / 310 *                | 
| IRON (mg)           | 26.81 | 15.13 |                8 / 18  *                 |
| SELENIUM (mcg)      | 124.85 | 111.17 |                 55                       |
| ZINC (mcg)          | 15.69 | 12.6 |                11 / 8  *                 |
| VITAMIN B12 (mcg)   | 6.34 | 4.97 |                 2.4                      |
| VITAMIN A (mcg)     | 797.19 | 951.51 |               900 / 700                  |
| VITAMIN C (mg)      | 129.21   | 89.64   |                90 / 75                   |
| VITAMIN D (mcg)     | 8.59 | 5.46 |                 15 *                     |
| VITAMIN E (mg)      | 7.71 | 10.72 |                 15                       |
| VITAMIN K (ug)      | 154.07   | 228.19   |               120 / 90                   |


Please note, values marked with an asterisk are the recommended dietary allowance or adequate intake values for adults between 18 and 30 years old. These values are the daily intake level that is sufficient to meet the needs of 97-98% of the population. If you are older than 30, your recommended dietary allowance may be higher or lower than what's shown here. You can visit https://ods.od.nih.gov/Health_Information/Dietary_Reference_Intakes.aspx for a complete list of dietary reference intakes by age and gender.

If your intake of a vitamin or mineral falls below the recommended level for your gender you should consider increasing your intake of foods that are a good source of that vitamin or mineral. For example, the study average intake of vitamin D falls below the recommended level of 15 mcg. So, many of the participants in our study could benefit by eating more vitamin D containing foods like fatty fish (such as salmon, tuna, and mackerel), fortified milk products, beef liver, egg yolks, and mushrooms. 


\newpage
##Figure 2:
![](/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-9-1.pdf)<!-- --> 

**Figure 2** shows your fiber intake each day throughout the study. You had an average fiber consumption of 23.53 grams per day. The red line in this plot shows the recommended level of fiber consumption in grams, while the gold line shows your average daily fiber intake. The current recomendation for fiber intake is to consume 14g per 1000kcal. We calculated your fiber needs to be 31 grams per day.

Fiber is an important energy source for the bacteria that live within your microbiome. Meeting the recommended intake of fiber may help to support a healthy microbiome. If your average fiber intake falls below the recommended intake level, you should consider increasing your intake of fiber by eating more fibrous vegetables and/or increasing your intake of whole grains.


##Figure 3:
![](/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-10-1.pdf)<!-- --> 

**Figure 3** shows the relative abundance of the major food groups you consumed on a day to day basis during the study, excluding water. Each color is representative of a unique food group as annotated on the figure legend. Please note that these food groups are displayed by their weight in grams, so the plot tends to over-represent the contribution of liquid foods and foods with high water content to your diet. You may notice that the largest food component each day is called "Sugars, sweets and beverages". For most participants in this study, the size of that bar is driven by intake of coffee and tea, not sugar and candy.


\newpage
#Your Microbiome Composition Data

The next series of plots show information about the composition of your microbiome.  

We sequenced the DNA of the microbes that live within your gut from the daily stool samples you provided. We used a type of sequencing called Shallow-shotgun sequencing coupled with programs developed in the Knights's lab to identify the species of bacteria that live inside you.


##Figure 4:
![](/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 

**Figure 4** shows the most abundant bacterial species within your gut per each day of the study. The "<2.5% abundance" column represents a sum of bacterial species that individually account for less than 2.5% of the bacterial population within your gut microbiome. If you are missing microbiome data for any day on this plot it's either becuase you did not submit a sample for that day, your sample wasn't able to be sequenced, or we did not get enough sequence data from your sample to reliably use it for data analysis.



##Figure 5:
![](/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-12-1.pdf)<!-- --> 

**Figure 5** shows how the bacterial diversity exhibited within your gut changes on a daily basis. The measure of diversity shown here is called the Shannon index of alpha diversity. The Shannon index accounts for both the abundance and eveness of bacterial species present within the gut microbiome. A higher Shannon index typically indicates a more diverse microbiome community.

\newpage
##Figure 6:

![](/Users/abby/Documents/Projects/zebra/output/output/report.MCTs01_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 
**Figure 6** is a PCoA plot. A PCoA plot summarizes variability within a given dataset by producing a set of uncorrelated axes. Essentially, the PCoA plot can be utilized to interpret similarity of data points -- data points closer to one another are more similar to one another, while points further away from eachother are more dissimilar. Each data point in this plot represents a unique sample collected one of the participants in this study. Samples from the same study participant are the same color. 

The points corresponding to your microbiome are outlined in black and are highlighted by the smaller ellipse. When viewing this plot, keep in mind that it is a 2-dimensional representation of a 3-dimensional plot. So if your points look like they are all in the same place, they might actually extend into or out of the page!


\newpage
#Your blood draw data

The results of your blood draw at the start of the dietary supplementation period are shown in the table below along with the reference ranges. If any of your values fall outside the reference ranges you should contact your health care provider to have these test re-run.


| Test                  | Result                              | Normal Range |
|-----------------------|-------------------------------------|--------------|
| Cholesterol (mg/dL)   | 165   |  <200        |
| Triglycerides (mg/dL) | 65         |  <150        | 
| HDL (mg/dL)           | 72           |  >60         |  
| LDL (mg/dL)           | 81           |  <100        |   
| Glucose (mg/dL)       | 87           | 70 - 100     |
