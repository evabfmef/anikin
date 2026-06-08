"""
# Machine Learning Classification of Kin Discrimination Groups

**Author:** Eva Stare

This script applies supervised machine learning to binary presence–absence matrices
(six genomic feature classes) to identify features predictive of kin discrimination
group membership.

**Note:** The same pipeline was applied independently to both the 39-strain
(14 KD groups) and 67-strain datasets; only the input file differs.

**Groups definition**
Groups that contained 2 or less strains were all assigned to KDoutgroups. The strains in the KDoutgroups category were not kin strains. The kin groups they belonged to had 2 or fewer representatives, leading us to combine them into the larger KDoutgroups category for our ML analysis. This strategy was chosen to build a more substantial dataset, thus enhancing the effectiveness of the ML process.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_validate, cross_val_score
from sklearn.metrics import make_scorer, accuracy_score, f1_score

from pprint import pprint

import importlib

for pkg in ["pandas", "numpy", "matplotlib", "seaborn", "sklearn", "joblib"]:
    try:
        mod = importlib.import_module(pkg)
        print(f"{pkg}: {mod.__version__}")
    except Exception as e:
        print(f"{pkg}: ERROR ({e})")

"""
## Load data
"""

file_path = '67_megamatrix_for_ML.xlsx'
data = pd.read_excel(file_path, index_col=0)

X = data.drop(columns=['label_gruped', 'lab_label'])
y = data['label_gruped']

"""
What it does:
- Loads your Excel file containing the feature presence/absence matrix
- index_col=0 means the first column becomes the row index (strain names)
- Creates X (features): All columns except the labels — this is your binary feature matrix (0/1 values)
- Creates y (target): The label_gruped column containing kin group labels (G1, G4, G9, G11, KDoutgroups)
"""

"""
## Define models and cross-validation
"""

scorers = {'accuracy': make_scorer(accuracy_score),
           'f1': make_scorer(f1_score, average='weighted')}


# Initialize classifiers with default parameters
rf_model = RandomForestClassifier()
svm_model = SVC()
nb_model = GaussianNB()
dummy_clf = DummyClassifier(strategy="most_frequent")


# Perform 3-fold cross-validation for all models
models = {
    'Random Forest': rf_model,
    'SVM': svm_model,
    'Naive Bayes': nb_model,
    'Dummy': dummy_clf
}


# Create a StratifiedKFold object to maintain the proportion of the target variable
stratified_kf = StratifiedKFold(n_splits=3)
cv_results = {}

for model_name, model in models.items():
    cv_results[model_name] = cross_validate(model, X, y, cv=stratified_kf, scoring=scorers)

"""
#### 1. Defines scoring metrics:
- **accuracy**: Proportion of correct predictions
- **f1** (weighted): Harmonic mean of precision and recall, weighted by class frequency — important for imbalanced classes


#### 2. Initializes four classifiers:
- **Random Forest**: Builds many decision trees and averages their predictions. Good for high-dimensional data, handles feature interactions.
- **SVM (Support Vector Machine)**: Finds the optimal hyperplane separating classes. Default uses RBF kernel.
- **Naive Bayes**: Assumes features are independent given the class. Fast and simple.
- **Dummy**: Always predicts the most frequent class — this is your baseline. If real models don't beat this, they're not learning anything useful.


#### 3. Sets up 3-fold stratified cross-validation:
- Data is split into 3 parts (folds)
- "Stratified" means each fold has the same proportion of kin groups as the full dataset
- Each model is trained on 2 folds and tested on the remaining 1, rotating through all combinations

#### 4. Runs cross-validation for each model and stores results
"""

"""
## Calculate and display results
"""

# Averaging F1 and Accuracy scores over all folds for each model
average_results = {}
for model_name, scores in cv_results.items():
    average_accuracy = scores['test_accuracy'].mean()
    average_f1 = scores['test_f1'].mean()
    average_results[model_name] = {
        'Average Accuracy': average_accuracy,
        'Average F1 Score': average_f1}

model_dfs = []
for model, scores in cv_results.items():
    avg_accuracy = average_results[model]['Average Accuracy']
    std_accuracy = np.std(scores['test_accuracy'])
    avg_f1 = average_results[model]['Average F1 Score']
    std_f1 = np.std(scores['test_f1'])

    model_df = pd.DataFrame({
        "Model": [model],
        "Accuracy": [f"{avg_accuracy:.3f} (+/- {std_accuracy:.3f})"],
        "F1 Score": [f"{avg_f1:.3f} (+/- {std_f1:.3f})"]
    })

    model_dfs.append(model_df)

# Concatenating all model DataFrames
formatted_average_results_concat = pd.concat(model_dfs).set_index("Model")
pprint(formatted_average_results_concat)

"""
For each model, it calculates:
- Mean accuracy across 3 folds
- Mean F1 score across 3 folds
- Standard deviation — indicates how consistent the model is across folds

Creates a formatted table showing mean (+/- std) for each metric.

If standard deviation is low, it suggests consistent performance across folds.

---
**Results**:

- **Random Forest** performed best, followed by Naive Bayes
- All real models vastly outperform the dummy baseline
"""

"""
## Visualization (boxplots)
"""

accuracy_scores_for_plot = {model: scores['test_accuracy'] for model, scores in cv_results.items()}
f1_scores_for_plot = {model: scores['test_f1'] for model, scores in cv_results.items()}

# Convert to DataFrame for seaborn compatibility
accuracy_scores_df = pd.DataFrame(accuracy_scores_for_plot)
f1_scores_df = pd.DataFrame(f1_scores_for_plot)

# Setting up the plot
plt.figure(figsize=(14, 6))

# Plot for Accuracy
plt.subplot(1, 2, 1)
sns.boxplot(data=accuracy_scores_df)
plt.title('Accuracy Scores')
plt.ylabel('Accuracy')
plt.xlabel('Model')

# Plot for F1 Score
plt.subplot(1, 2, 2)
sns.boxplot(data=f1_scores_df)
plt.title('F1 Scores')
plt.ylabel('F1 Score')
plt.xlabel('Model')

plt.tight_layout()
plt.show()

"""
What it does:

- Extracts accuracy and F1 scores from each fold for each model
- Creates side-by-side boxplots showing the distribution of scores

Why boxplots? They show:

- The median (line in the box)
- Spread of scores (box height)
- Any outliers
- Visual comparison between models

With only 3 data points per model (one per fold), the boxplots are simple, but they still show relative performance and consistency.
"""

"""
## Save all models
"""

import joblib

# Train and save all models on full data (excluding Dummy)
models_to_save = {
    'Random Forest': RandomForestClassifier(),
    'Naive Bayes': GaussianNB(),
    'SVM': SVC()
}

for model_name, model in models_to_save.items():
    model.fit(X, y)
    safe_name = model_name.lower().replace(' ', '_')
    save_path = f'KD_model_{safe_name}.joblib'
    joblib.dump(model, save_path)
    print(f"Saved: {save_path}")

"""
This chunk:
- Creates fresh model instances and trains each on ALL data
- Saves models to disk using joblib:
  - KD_model_random_forest.joblib
  - KD_model_naive_bayes.joblib
  - KD_model_svm.joblib

"""

n_features = X.shape[1]
print(f"Total number of features: {n_features}")

"""
---
# Feature Importance Analysis

We use **Naive Bayes** and **Random Forest** — our two best-performing models — to identify features that discriminate between kin groups.
"""

"""
## Naive Bayes discriminative analysis
"""

nb_model = GaussianNB()
nb_model.fit(X, y)

class_means = nb_model.theta_

# Create class_means_df (needed for pattern analysis later)
class_means_df = pd.DataFrame(class_means, columns=X.columns, index=nb_model.classes_)

# Use variance to measure discriminative power
mean_variance = np.var(class_means, axis=0)

nb_features_df = pd.DataFrame({
    'Feature': X.columns,
    'Discriminative_Power': mean_variance
})
nb_features_df_sorted = nb_features_df.sort_values(by='Discriminative_Power', ascending=False)
nb_features_df_sorted.head(20)

len(nb_model.classes_)

"""
What it does:

1. Trains a Naive Bayes model
2. Extracts *theta_* — the mean of each feature for each class
    - Shape: (7 kin groups × ~4800 features)
    - For binary data, this is essentially the proportion of strains in each group that have the feature

3. Calculates variance across classes for each feature
    - High variance = feature presence differs a lot between kin groups = discriminative
    - Low variance = feature is similarly present/absent across all groups = not useful for classification

4. Sorts features by discriminative power

Why variance instead of range?
With 7 groups and binary data:
- Range (max - min) maxes out at 1.0 for many features
- Variance captures more nuance
"""

"""
**Result:**
The 0.240009 is the maximum possible variance for your data.
Because we have 7 kin groups. Variance is maximized when data is split evenly (3/7 or 4/7 split - either 3 groups with the feature, 4 without, or vice versa — both give the same variance of 0.24).
"""

"""
## View class patterns
"""

# Check your classes
print("Kin groups:", nb_model.classes_)
print("Number of groups:", len(nb_model.classes_))

# Look at the actual mean pattern for top features
top_features = nb_features_df_sorted.head(20)['Feature'].tolist()
class_means_df[top_features]

"""
What it does:

1. Shows you which kin groups exist: G25_A' 'G50' 'G53' 'G60_A' 'G60_B' 'G63_A' 'KDoutgroups'

2. Creates a pattern table showing, for each top feature:
- Row = kin group
- Value = proportion of strains in that group with the feature (0.0 to 1.0)

Example:
```
                    feature_X    feature_Y
G1                    1.0       0.0
G11                   1.0       1.0
G4                    0.0       1.0
G9                    0.0       1.0
KDoutgroups           0.0       0.92
```
feature_X: Present in G1 and G11, absent elsewhere — could be a G1/G11 marker    
feature_Y: Absent only in G1 — could help identify non-G1 strains
"""

"""
## Export NB pattern table
"""

# Export the full pattern table for top 100 features of NB
top_100_features = nb_features_df_sorted.head(100)['Feature'].tolist()
pattern_table = class_means_df[top_100_features].T  # transpose so features are rows
pattern_table['Discriminative_Power'] = nb_features_df_sorted.head(100)['Discriminative_Power'].values
pattern_table.to_excel('feature_patterns_by_kin_group_NB.xlsx')
print("Exported: feature_patterns_by_kin_group_NB.xlsx")

"""
What it does:

1. Takes top 100 discriminative features
2. Creates a table with:
    - Rows = features
    - Columns = kin groups + discriminative power score

3. Exports to Excel

We get a spreadsheet where we can see which features define which groups and sort/filter by group.
"""

"""
## Random Forest feature importances
"""

# Feature importances (only from Random Forest)
rf_model = models_to_save['Random Forest']
feature_importances = rf_model.feature_importances_

features_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': feature_importances
})
features_df_sorted = features_df.sort_values(by='Importance', ascending=False)
features_df_sorted.head(100).to_excel('feature_importances.xlsx', index=False)
features_df_sorted[:100]

n_features = X.shape[1]
print(f"Total number of features: {n_features}")

"""
Extracts the *feature_importances_* attribute from the trained Random Forest model. This is an array of 4779 values (one per feature).

**How RF calculates importance:** For each feature, it measures how much that feature helped reduce "impurity" (classification error) across all trees in the forest. Features that create clean splits (one class on each side) get higher scores.

Random Forest importance values sum to 1.0 across all features. With 4779 features, if importance were distributed equally, each feature would have: 1.0 / 4779 ≈ **0.000209**.

Plus:
- Creates a DataFrame pairing each gefeature ne name with its importance score.
- Sorts features from most to least important.
- Exports top 100 features to Excel.
- Displays them in the notebook.
"""

"""
**Result:**

The Random Forest model identified **hypothetical protein_group_8849** as the top feature, with an importance of 0.0164, corresponding to ~1.64% of the total predictive power. Under a null model of equal contribution across all 4779 features, each feature would be expected to contribute 1/4779 ≈ 0.000209. Thus, the top feature alone contributes ~78× more signal than expected by chance, making it the most informative predictor in the model.

Our top 5 features (importance 0.0117 to 0.0164) together account for 0.068, or about **6.78%** of the total predictive power. Given the large feature space, this represents a substantial concentration of signal: under equal importance, the top five features would be expected to contribute only 5/4779 ≈ 0.00105 (≈ 0.10%) of the total importance. Therefore, the observed contribution of the top five features reflects an enrichment of roughly ~65× over random expectation.

For genomic data with thousands of features, these are reasonable values.
Feature importance values were interpreted relative to the null expectation of equal contribution across all features (1/4779 ≈ 0.000209). Importance values exceeding this baseline by one or more orders of magnitude were considered enriched and indicative of non-random predictive contribution.

The highest-ranked features exhibited importance values ~78× greater than expected under equal importance, indicating a strong concentration of predictive signal despite the high dimensionality of the feature space.

Importance Interpretation
- ≥ 0.01 Strong contributor (your top features)
- 0.005 0.01 → Moderate contributor
- 0.001 0.005 → Weak but real signal
< 0.001 Likely noise
"""

print(f"Top 10 account for: {features_df_sorted.head(10)['Importance'].sum():.1%}")
print(f"Top 50 account for: {features_df_sorted.head(50)['Importance'].sum():.1%}")
print(f"Top 100 account for: {features_df_sorted.head(100)['Importance'].sum():.1%}")

"""
### Compare top genes from RF and NB
"""

# Check the top 10 from each model
print("RF top 10:")
print(features_df_sorted.head(10)['Feature'].tolist())

print("\nNB top 10:")
print(nb_features_df_sorted.head(10)['Feature'].tolist())

"""
**Results:**

**RF top 10**:
- hypothetical protein_group_8849
- **hypothetical protein_ydbD** (general stress protein, survival of ethanol stress)
- **yrdK** (Unknown, regulated by: AzlB regulon)
- hypothetical protein_group_9127
- **yydC.0** (Unknown, regulated by: SigB regulon)
- **yqgA** (unknown, genetically related to the cell wall-degrading dl-endopeptidases)
- hypothetical protein_group_9702
- **desA, des, yocE** (phospholipid desaturase; adaptation of membrane fluidity at low temperatures)
- **rpmGA** (ribosomal protein bL33a; non-essential; translation)
- **yxnB** (Unknown)

**NB top 10** (**NB top 10** (see above for all features including the highest discriminative power!)):
- hypothetical protein_group_4237
- hypothetical protein_group_9692
- yopF.1
- hypothetical protein_group_11347
- hypothetical protein_group_659
- hypothetical protein_group_11348
- hypothetical protein_group_11346
- hypothetical protein_group_11345
- hypothetical protein_group_11344
- hypothetical protein_group_11343


Key insight: Completely different genes! The methods are capturing different aspects of what makes genes useful for classification.

**How Naive Bayes picks genes**

NB looks at each gene **independently** and asks: *"How different is this gene's presence across kin groups?"*

Example:
```
           Gene_A
G1           1.0
G11          1.0
G4           0.0
G9           0.0
KDoutgroups  0.0
```
Present in G1 and G11, absent elsewhere. High discriminative power.

**How Random Forest picks genes**

RF builds decision trees and asks: *"Which gene creates the best split at this point in the tree, given what other genes already did?"*

Example:
```
           Gene_A    Gene_B    Gene_C
G1           1.0       1.0       0.0
G11          1.0       0.0       1.0
G4           0.0       1.0       0.0
G9           0.0       0.0       1.0
KDoutgroups  0.0       0.5       0.5
```

Gene_A alone separates G1+G11 from the rest. But to separate G1 from G11, you need Gene_B or Gene_C.

RF might rank Gene_B and Gene_C highly because they're useful in combination with Gene_A. But NB would rank them lower because individually their patterns aren't as clean.

If you have 50 genes with nearly identical patterns (for example all part of the same prophage region):
```
           Gene_1   Gene_2   Gene_3  ... Gene_50
G1           1.0      1.0      1.0        1.0
G11          1.0      1.0      1.0        1.0
G4           0.0      0.0      0.0        0.0
G9           0.0      0.0      0.0        0.0
KDoutgroups  0.0      0.0      0.0        0.0
```

**Naive Bayes**: All 50 genes get high discriminative power (0.24 each) — they all have the same clean pattern.
**Random Forest**: Only a few of these 50 get high importance. Because once the tree uses Gene_1 to make a split, Gene_2-50 don't add new information. RF spreads the importance thinly across redundant genes, or picks just one.

Summary

```
Method           What it rewards                         Weakness
Naive Bayes      Clean individual patterns               Ranks redundant genes equally high
Random Forest    Unique contribution to classification   Dilutes importance across similar genes
``
"""

"""
## Find overlap between methods
"""

# Compare top 100 from each method
rf_top100 = set(features_df_sorted.head(100)['Feature'])
nb_top100 = set(nb_features_df_sorted.head(100)['Feature'])

print(f"Overlap: {len(rf_top100 & nb_top100)} genes")
print(f"Shared genes: {rf_top100 & nb_top100}")

"""
What it does: Finds genes that appear in the top 100 for BOTH methods.

Results: Only **14** genes overlap!

These genes passed both filters:

Clean presence/absence patterns (NB scores high)
Unique predictive contribution (RF scores high)
Why so little overlap?

Random Forest finds genes useful in combination with others (interaction effects)
Naive Bayes finds genes with the cleanest individual presence/absence patterns
Different approaches to the same classification problem


**Genes:**
- **group_7478**
- hypothetical protein_group_7476
- **nisC** (nisin biosynthesis)
- **yrkO** (Unknown membrane protein, regulated by YrkP regulon)
- hypothetical protein_group_1853
- hypothetical protein_group_2516
- **group_5168**
- hypothetical protein_group_659
- hypothetical protein_group_9700
- hypothetical protein_group_4237
- hypothetical protein_group_11345
- **yopF.1** (Unknown, SP-beta prophage product)
- hypothetical protein_group_7475
- hypothetical protein_group_7650
"""

"""
## Rank the overlapping genes
"""

overlap_genes = rf_top100 & nb_top100

for gene in overlap_genes:
    rf_rank = (features_df_sorted['Feature'] == gene).argmax() + 1
    nb_rank = (nb_features_df_sorted['Feature'] == gene).argmax() + 1
    print(f"{gene}: RF rank {rf_rank}, NB rank {nb_rank}")

"""
### If you want genes that cleanly define each kin group → Naive Bayes

NB directly shows you which features have distinct presence/absence patterns per group. The pattern table tells you exactly:
- "This feature is present in G1+G11, absent elsewhere"
- "This feature is only in G9"

This is intuitive and directly interpretable for biology. You can look at a feature and immediately understand which groups it distinguishes.

### If you want the minimal set of features needed to classify strains → Random Forest

RF tells you which features are actually needed for accurate prediction. If 50 features have identical patterns, RF will (roughly) pick one or spread importance across them, while NB ranks all 50 equally high.



---


**We can use both:**

1. Start with Naive Bayes to identify candidate features with clean patterns — these are easier to interpret. The pattern table (gene_patterns_by_kin_group_NB.xlsx) is directly useful for biological interpretation.
2. Then check RF importance to filter out redundancy. If NB gives you 100 top features but many have identical patterns (likely from the same genomic region), RF can help you identify which ones are actually contributing unique information.
3. Prioritize the overlapping features — both methods agree these matter, making them your strongest candidates.
"""

"""
## Find robust feature candidates

Now it takes the top 200 features from each method and converts them to sets (unordered collections of unique items) cause with top 100, we only got 14 overlapping features.
Expanding to 200 catches more features that both methods agree are important, giving you 39 candidates.
"""

# Find genes NB says are discriminative AND RF says are important
# Use top 200 from each to catch more overlap
rf_top200 = set(features_df_sorted.head(200)['Feature'])
nb_top200 = set(nb_features_df_sorted.head(200)['Feature'])

robust_genes = rf_top200 & nb_top200
print(f"Robust candidates: {len(robust_genes)} genes")

# Get their patterns
robust_list = list(robust_genes)
class_means_df[robust_list].round(2)

"""
What it does:
- The & operator finds the intersection — features that appear in BOTH sets. These are your robust candidates because two completely different methods independently identified them as important.
- Converts the set to a list (needed for indexing the DataFrame), then displays the presence patterns for these features across all kin groups, rounded to 2 decimal places.
"""

overlap_genes = rf_top100 & nb_top100

for gene in overlap_genes:
    rf_rank = (features_df_sorted['Feature'] == gene).argmax() + 1
    nb_rank = (nb_features_df_sorted['Feature'] == gene).argmax() + 1
    print(f"{gene}: RF rank {rf_rank}, NB rank {nb_rank}")

# Export the robust candidates with their patterns and rankings
robust_list = list(robust_genes)

# Get patterns
robust_patterns = class_means_df[robust_list].T

# Add RF and NB rankings
robust_patterns['RF_rank'] = [
    (features_df_sorted['Feature'] == gene).argmax() + 1
    for gene in robust_list
]
robust_patterns['NB_rank'] = [
    (nb_features_df_sorted['Feature'] == gene).argmax() + 1
    for gene in robust_list
]

# Sort by average rank
robust_patterns['Avg_rank'] = (robust_patterns['RF_rank'] + robust_patterns['NB_rank']) / 2
robust_patterns = robust_patterns.sort_values('Avg_rank')

robust_patterns.to_excel('robust_candidates.xlsx')
robust_patterns

"""
### Results: Genomic features defining kin discrimination groups

Using machine-learning–based feature selection on a presence/absence matrix of genomic features across Bacillus subtilis strains assigned to kin discrimination (KD) groups, we identified a set of strain-variable loci that consistently contributed to KD group separation. Although hypothetical proteins remain a substantial component of the signal, the updated feature set reveals clearer enrichment for **cell envelope biosynthesis**, **mobile genetic elements**, **regulatory systems**, and **developmental pathways**, highlighting specific biological mechanisms underlying kin differentiation.

Overall, the selected features are strongly biased toward non-core, horizontally variable, and adaptive genomic components, rather than conserved housekeeping genes. This pattern indicates that kin discrimination groups are structured primarily by flexible genomic regions associated with surface properties, horizontal gene transfer, and ecological interactions.

Major functional themes
1. **Persistent dominance of hypothetical and unknown proteins**

A large fraction of KD-defining features lack detailed functional annotation. Multiple hypothetical protein groups (e.g. hypothetical protein_group_2085, _843, _2140, _2598, _73, _2154, _1869, _2734) and several completely uncharacterized genes (e.g. yxzG, yyaO, yrkN, yxiF, yobF) appear among the most informative predictors.

This continued enrichment suggests that poorly annotated, strain-specific proteins represent a major reservoir of kin-group–specific functions, likely including novel membrane proteins, toxins, or regulatory elements that have not yet been functionally characterized.

2. **Strong signal from cell envelope and teichoic acid biosynthesis**

A newly prominent theme is the enrichment of genes involved in teichoic acid and lipoteichoic acid synthesis, including:

- tagA variants (teichoic acid biosynthesis)
- tarQ-like glucosyltransferase (cell wall modification)
- multiple lipoteichoic acid acylation reactions (MNXR136290, MNXR136296, MNXR136300, MNXR136302)

Teichoic acids are major determinants of cell surface charge, structure, and interaction specificity, strongly suggesting that cell-surface chemistry is a key axis of kin discrimination. Variation in these pathways likely alters recognition, adhesion, susceptibility to antagonistic systems, and contact-dependent interactions.

3. **Continued enrichment of prophage and mobile element features**

Prophage-associated genes are important KD discriminators, including:

- SP-β prophage membrane proteins (yopD, yopE)
- prophage encoded toxin-antixoin system (BF11) (PBSX-associated toxin–antitoxin system

This persistent signal supports a model in which horizontal gene transfer and prophage cargo genes shape social group structure, potentially by introducing toxins, immunity functions, or regulatory modules that drive compatibility between strains.

4. **Competition and immunity-related systems**

Several features are directly linked to inter-strain antagonism and protection:

- wapI (immunity against the toxic activity of WapA)
- prophage encoded toxin-antixoin system (BF11) (PBSX-associated toxin–antitoxin system

Membrane-associated prophage proteins (yopD, yopE) are also likely involved in host interaction or defense

This reinforces the idea that kin discrimination is tightly coupled to competitive interactions, where shared immunity and antagonistic systems define social compatibility.

5. **Regulatory and signaling components**

Regulatory elements also contribute to KD group structure, including:

- phrK peptide involved in ComA-dependent quorum sensing and competence regulation
- LysR-family transcriptional regulator yxjO

These results indicate that signaling pathways and transcriptional rewiring accompany genomic differentiation, likely coordinating prophage activity, surface remodeling, and stress responses in a group-specific manner.

6. **Developmental and sporulation-linked genes**

Multiple newly identified sporulation-associated genes (yyaL, ykzH) are enriched among KD-defining features, suggesting that life-history strategies and developmental timing differ systematically between kin groups. Such variation could influence biofilm formation, persistence, and competitive outcomes.

7. **Continued depletion of core housekeeping functions**

As in the earlier analysis, very few conserved metabolic or essential genes appear among KD-defining loci. The signal remains dominated by accessory and flexible genome components, emphasizing that kin discrimination is encoded primarily in the variable genome.

**Summary**

Together, the updated feature set confirms and refines the earlier conclusions: kin discrimination groups in Bacillus subtilis are primarily defined by adaptive, horizontally variable genomic components, rather than by conserved core biology. However, the new results reveal a particularly strong mechanistic emphasis on:

- Cell envelope architecture (teichoic and lipoteichoic acid biosynthesis)

- Mobile genetic elements and prophage cargo

- Competition and immunity systems

- Regulatory and developmental pathways

This pattern strongly supports a model in which kin recognition and social compatibility emerge from rapidly evolving surface properties, mobile elements,
and antagonistic systems, enabling fine-scale ecological differentiation among closely related strains.
"""
