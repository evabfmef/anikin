"""
Machine Learning Classification of Kin Discrimination Groups

Author: Eva Stare

Applies supervised machine learning to binary presence/absence matrices (six
genomic feature classes) to identify features predictive of kin discrimination
(KD) group membership.

The same pipeline was applied independently to the 39-strain (14 KD groups)
and 67-strain datasets; only the input file differs.

Group definition
----------------
KD groups containing <=2 strains were pooled into a single "KDoutgroups"
category to build a more substantial dataset for the ML analysis. Strains in
KDoutgroups are not kin to one another.

Outputs
-------
  KD_model_random_forest.joblib, KD_model_naive_bayes.joblib, KD_model_svm.joblib
  feature_patterns_by_kin_group_NB.xlsx
  feature_importances.xlsx
  robust_candidates.xlsx

Environment
-----------
See requirements.txt for pinned package versions.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import importlib
from pprint import pprint

from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import make_scorer, accuracy_score, f1_score

# Report package versions for reproducibility
for pkg in ["pandas", "numpy", "matplotlib", "seaborn", "sklearn", "joblib"]:
    try:
        mod = importlib.import_module(pkg)
        print(f"{pkg}: {mod.__version__}")
    except Exception as e:
        print(f"{pkg}: ERROR ({e})")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# Set to the 39- or 67-strain feature matrix as appropriate.
file_path = '67_megamatrix_for_ML.xlsx'
data = pd.read_excel(file_path, index_col=0)   # first column = strain names

X = data.drop(columns=['label_gruped', 'lab_label'])  # binary feature matrix
y = data['label_gruped']                              # KD group labels


# ---------------------------------------------------------------------------
# Define models and cross-validation
# ---------------------------------------------------------------------------
scorers = {'accuracy': make_scorer(accuracy_score),
           'f1': make_scorer(f1_score, average='weighted')}

models = {
    'Random Forest': RandomForestClassifier(),
    'SVM': SVC(),
    'Naive Bayes': GaussianNB(),
    'Dummy': DummyClassifier(strategy="most_frequent"),  # most-frequent baseline
}

# Stratified 3-fold CV preserves class proportions in each fold.
stratified_kf = StratifiedKFold(n_splits=3)
cv_results = {}
for model_name, model in models.items():
    cv_results[model_name] = cross_validate(model, X, y, cv=stratified_kf, scoring=scorers)


# ---------------------------------------------------------------------------
# Calculate and display results
# ---------------------------------------------------------------------------
average_results = {}
for model_name, scores in cv_results.items():
    average_results[model_name] = {
        'Average Accuracy': scores['test_accuracy'].mean(),
        'Average F1 Score': scores['test_f1'].mean(),
    }

model_dfs = []
for model, scores in cv_results.items():
    avg_accuracy = average_results[model]['Average Accuracy']
    std_accuracy = np.std(scores['test_accuracy'])
    avg_f1 = average_results[model]['Average F1 Score']
    std_f1 = np.std(scores['test_f1'])
    model_dfs.append(pd.DataFrame({
        "Model": [model],
        "Accuracy": [f"{avg_accuracy:.3f} (+/- {std_accuracy:.3f})"],
        "F1 Score": [f"{avg_f1:.3f} (+/- {std_f1:.3f})"],
    }))

formatted_average_results_concat = pd.concat(model_dfs).set_index("Model")
pprint(formatted_average_results_concat)


# ---------------------------------------------------------------------------
# Visualisation (boxplots of per-fold scores)
# ---------------------------------------------------------------------------
accuracy_scores_df = pd.DataFrame(
    {model: scores['test_accuracy'] for model, scores in cv_results.items()})
f1_scores_df = pd.DataFrame(
    {model: scores['test_f1'] for model, scores in cv_results.items()})

plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
sns.boxplot(data=accuracy_scores_df)
plt.title('Accuracy Scores')
plt.ylabel('Accuracy')
plt.xlabel('Model')

plt.subplot(1, 2, 2)
sns.boxplot(data=f1_scores_df)
plt.title('F1 Scores')
plt.ylabel('F1 Score')
plt.xlabel('Model')

plt.tight_layout()
plt.savefig('ml_model_comparison.png', dpi=150, bbox_inches='tight')
plt.show()


# ---------------------------------------------------------------------------
# Train and save final models on the full dataset (Dummy excluded)
# ---------------------------------------------------------------------------
models_to_save = {
    'Random Forest': RandomForestClassifier(),
    'Naive Bayes': GaussianNB(),
    'SVM': SVC(),
}
for model_name, model in models_to_save.items():
    model.fit(X, y)
    safe_name = model_name.lower().replace(' ', '_')
    save_path = f'KD_model_{safe_name}.joblib'
    joblib.dump(model, save_path)
    print(f"Saved: {save_path}")

n_features = X.shape[1]
print(f"Total number of features: {n_features}")


# ===========================================================================
# Feature importance analysis (Naive Bayes + Random Forest)
# ===========================================================================

# ---------------------------------------------------------------------------
# Naive Bayes discriminative analysis
# ---------------------------------------------------------------------------
# theta_ is the per-class feature mean (for binary data, the proportion of
# strains in each group carrying the feature). Variance of theta_ across
# classes is used as a discriminative-power score.
nb_model = GaussianNB()
nb_model.fit(X, y)

class_means = nb_model.theta_
class_means_df = pd.DataFrame(class_means, columns=X.columns, index=nb_model.classes_)

mean_variance = np.var(class_means, axis=0)
nb_features_df = pd.DataFrame({
    'Feature': X.columns,
    'Discriminative_Power': mean_variance,
})
nb_features_df_sorted = nb_features_df.sort_values(by='Discriminative_Power', ascending=False)

print("Kin groups:", nb_model.classes_)
print("Number of groups:", len(nb_model.classes_))

# Export NB pattern table for the top 100 discriminative features
top_100_features = nb_features_df_sorted.head(100)['Feature'].tolist()
pattern_table = class_means_df[top_100_features].T   # features as rows
pattern_table['Discriminative_Power'] = nb_features_df_sorted.head(100)['Discriminative_Power'].values
pattern_table.to_excel('feature_patterns_by_kin_group_NB.xlsx')
print("Exported: feature_patterns_by_kin_group_NB.xlsx")


# ---------------------------------------------------------------------------
# Random Forest feature importances
# ---------------------------------------------------------------------------
rf_model = models_to_save['Random Forest']
feature_importances = rf_model.feature_importances_

features_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': feature_importances,
})
features_df_sorted = features_df.sort_values(by='Importance', ascending=False)
features_df_sorted.head(100).to_excel('feature_importances.xlsx', index=False)
print("Exported: feature_importances.xlsx")

print(f"Top 10 account for:  {features_df_sorted.head(10)['Importance'].sum():.1%}")
print(f"Top 50 account for:  {features_df_sorted.head(50)['Importance'].sum():.1%}")
print(f"Top 100 account for: {features_df_sorted.head(100)['Importance'].sum():.1%}")


# ---------------------------------------------------------------------------
# Compare top features between methods
# ---------------------------------------------------------------------------
print("\nRF top 10:")
print(features_df_sorted.head(10)['Feature'].tolist())
print("\nNB top 10:")
print(nb_features_df_sorted.head(10)['Feature'].tolist())

# Overlap of the top 100 from each method
rf_top100 = set(features_df_sorted.head(100)['Feature'])
nb_top100 = set(nb_features_df_sorted.head(100)['Feature'])
overlap_genes = rf_top100 & nb_top100
print(f"\nTop-100 overlap: {len(overlap_genes)} genes")
print(f"Shared genes: {overlap_genes}")

for gene in overlap_genes:
    rf_rank = (features_df_sorted['Feature'] == gene).argmax() + 1
    nb_rank = (nb_features_df_sorted['Feature'] == gene).argmax() + 1
    print(f"  {gene}: RF rank {rf_rank}, NB rank {nb_rank}")


# ---------------------------------------------------------------------------
# Robust feature candidates (intersection of top 200 from each method)
# ---------------------------------------------------------------------------
# Top-100 overlap was small, so the intersection is taken over the top 200
# of each method to recover more mutually supported candidates.
rf_top200 = set(features_df_sorted.head(200)['Feature'])
nb_top200 = set(nb_features_df_sorted.head(200)['Feature'])
robust_genes = rf_top200 & nb_top200
print(f"\nRobust candidates (top-200 overlap): {len(robust_genes)} genes")

robust_list = list(robust_genes)

# Presence patterns for the robust candidates, annotated with RF/NB ranks.
robust_patterns = class_means_df[robust_list].T
robust_patterns['RF_rank'] = [
    (features_df_sorted['Feature'] == gene).argmax() + 1 for gene in robust_list]
robust_patterns['NB_rank'] = [
    (nb_features_df_sorted['Feature'] == gene).argmax() + 1 for gene in robust_list]
robust_patterns['Avg_rank'] = (robust_patterns['RF_rank'] + robust_patterns['NB_rank']) / 2
robust_patterns = robust_patterns.sort_values('Avg_rank')

robust_patterns.to_excel('robust_candidates.xlsx')
print("Exported: robust_candidates.xlsx")
