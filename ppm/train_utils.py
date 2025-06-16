""" Common functions used in training both canonical/cryptic and spliced models.
"""
import pandas as pd
from shap import TreeExplainer
from sklearn.model_selection import RandomizedSearchCV, LeaveOneGroupOut
import xgboost as xgb

PARAM_SETS = {
    "learning_rate"    : [0.1, 0.15, 0.20, 0.30 ] ,
    "max_depth"        : [ 2, 3, 4, 5, 6, 7],
    "min_child_weight" : [ 1, 2, 3],
    "gamma"            : [ 0.0, 0.1, 0.2, 0.3 ],
    "colsample_bytree" : [ 0.5, 0.7, 0.9, 0.95 ],
}

def add_shap_values(test_df, clf, feature_set):
    """ Function to add shap values describing the impact of a feature on the model outputs.
    """
    explainer = TreeExplainer(clf)
    shap_values = explainer.shap_values(test_df[feature_set])
    shap_df = pd.DataFrame(shap_values, columns=[f'{col}_shap' for col in feature_set])
    shap_df = shap_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)
    test_df = pd.merge(
        test_df,
        shap_df, how='inner', left_index=True, right_index=True
    )
    return test_df

def run_hpt(train_df, feature_set):
    """ Function to run hyperparameter tuning jobs.
    """
    clf = xgb.XGBClassifier()
    cv = RandomizedSearchCV(
        clf,
        param_distributions=PARAM_SETS, n_iter=10,
        scoring='balanced_accuracy', n_jobs=16,
        cv=LeaveOneGroupOut(), random_state=42, verbose=0,
        return_train_score=True,
    )
    searched_cv = cv.fit(train_df[feature_set], train_df['label'], groups=train_df['cvGroup'])
    return searched_cv.best_params_
