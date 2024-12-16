import pandas as pd
import warnings
# Suppress specific FutureWarnings from xgboost
warnings.filterwarnings('ignore', category=FutureWarning, 
                        message=".*is_sparse is deprecated.*")
warnings.filterwarnings('ignore', category=FutureWarning, 
                        message=".*is_categorical_dtype is deprecated.*")
import xgboost as xgb
from sklearn.model_selection import GridSearchCV, train_test_split, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression
import shap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from tqdm import tqdm

def prepare_data(data, target_col, covariates, seed = 42):
    data = data.dropna(subset=[target_col] + covariates, axis=0)
    # Extract all cases
    cases = data[data[target_col] == 1]
    
    # Randomly sample an equal number of controls
    controls = data[data[target_col] == 0].sample(n=len(cases), random_state=seed)
    
    # Concatenate cases and controls to form a balanced dataset
    balanced_data = pd.concat([cases, controls], axis=0).reset_index(drop=True)
    ## Print the distribution for the multimorbidity
    print(balanced_data.multimorbidity.value_counts())
    # Split the data into training (80%) and testing (20%) sets
    train_data, test_data = train_test_split(balanced_data, test_size=0.2, stratify=balanced_data[target_col], random_state=seed)
    covariates = [
        "GenMetS", "TVhour", "Sleephour", "PhysicalActivity",
        "SocialEconomicalStatus", "Age", "Education"
    ]    
    # Extract the features and target for each set
    X_train, y_train = train_data[covariates], train_data[target_col]
    X_test, y_test = test_data[covariates], test_data[target_col]
    X_train.reset_index(drop=True, inplace=True), y_train.reset_index(drop=True, inplace=True), X_test.reset_index(drop=True, inplace=True), y_test.reset_index(drop=True, inplace=True) 
    scaler = StandardScaler()
    trans_covariates  = [
    "GenMetS", "TVhour", "Sleephour", "PhysicalActivity",
    "SocialEconomicalStatus", "Age"
]
    X_train[trans_covariates] = scaler.fit_transform(X_train[trans_covariates])   
    X_test[trans_covariates] = scaler.transform(X_test[trans_covariates])
    return X_train, y_train, X_test, y_test
# Define the target column and covariates
target_col = 'multimorbidity'
covariates = [
    "GenMetS", "TVhour", "Sleephour", "PhysicalActivity",
    "SocialEconomicalStatus", "Age", "Education"
]

missing_covariates = ['GenMetS','Education']
# Load the data
# study_data_v2 = pd.read_csv('./study_data_v2.csv', sep = '\t')
study_data_v2 = pd.read_csv('./UKB_Asianwomenmen_s8792_mets_diseases_others.txt', sep='\t')
study_data_v2 = study_data_v2[study_data_v2.sex == 'Female']
## set eduction = NA when -3, -7
study_data_v2.loc[study_data_v2.Education == -3, 'Education'] = np.nan
study_data_v2.loc[study_data_v2.Education == -7, 'Education'] = np.nan

X_train, y_train, X_test, y_test = prepare_data(study_data_v2, target_col, covariates, seed=3)
def xgb_best(X_train, y_train):
    print("tunning for avoid overfitting")
    xgb1 =  xgb.XGBClassifier()
    parameters = {'objective':['binary:logistic'],
              'learning_rate': [0.01, 0.05, 0.1],       #default 0.1, lower better for avoid overfitting
              'max_depth': [2, 3, 4 ],               #default 3, lower better
              'min_child_weight': [1, 2, 4],        #default 1, larger better
              'subsample': [0.7,  1],           #default 1.0, lower better
              'colsample_bytree': [0.1, 0.5, 1.0],  #default 1.0, lower better
              }
    kfold = KFold(n_splits=5, shuffle=True)
    xgb1 = xgb.XGBClassifier(verbosity = 0, silent=True)
    parameters = {
        'objective':['binary:logistic'],
        'learning_rate': [0.01, 0.05, 0.1],  # default 0.1, lower better for avoid overfitting
        'max_depth': [2, 3, 4],  # default 3, lower better
        'min_child_weight': [1, 2, 4],  # default 1, larger better
        'subsample': [0.7, 1],  # default 1.0, lower better
        'colsample_bytree': [0.1, 0.5, 1.0],  # default 1.0, lower better
    }
    kfold = KFold(n_splits=5, shuffle=True)
    xgb_grid = GridSearchCV(xgb1,
                            parameters,
                            cv=kfold,
                            n_jobs=10,
                            verbose=True)
    xgb_grid.fit(X_train, y_train)
    best_pars = xgb_grid.best_params_

    xgb_best = xgb.XGBClassifier(
        colsample_bytree=best_pars['colsample_bytree'],
        learning_rate=best_pars['learning_rate'],
        max_depth=best_pars['max_depth'],
        gamma=5,
        min_child_weight=best_pars['min_child_weight'],
        subsample=best_pars['subsample'],
        verbosity = 0,
        silent=True
    )
    return xgb_best

def c_index(y_true, scores):
    n = len(y_true)
    assert len(scores) == n, "Input arrays should have the same length"
    
    concordant = 0
    permissible = 0
    tied = 0
    
    for i in range(n):
        for j in range(i+1, n):  # Compare all pairs of data points
            if y_true[i] != y_true[j]:  # Only consider pairs with differing outcomes
                permissible += 1
                if scores[i] < scores[j] and y_true[i] < y_true[j]:
                    concordant += 1
                elif scores[i] > scores[j] and y_true[i] > y_true[j]:
                    concordant += 1
                elif scores[i] == scores[j]:
                    tied += 1
    
    c_index = (concordant + 0.5*tied) / permissible
    return c_index


def train_evaluate_model(X_train, y_train, X_test, y_test, covariates, input_config = 'GenMetS', encoding_method='normalize'):
    if encoding_method not in ['dummy', 'normalize']:
        raise ValueError(f'Unknown encoding_method: {encoding_method}')

    if encoding_method not in ['dummy', 'normalize']:
        raise ValueError(f'Unknown encoding_method: {encoding_method}')

    # Adjusting feature selection based on input_config
    if input_config == 'GenMetS':
        X_train = X_train[['GenMetS']]
        X_test = X_test[['GenMetS']]
    elif input_config == 'cov + GenMetS':
        X_train = X_train[covariates]
        X_test = X_test[covariates]
    elif input_config == 'GenMetS + Age':
        # Include only 'GenMetS' and 'Age'
        X_train = X_train[['GenMetS', 'Age']]
        X_test = X_test[['GenMetS', 'Age']]
    elif input_config == 'cov':
        # Include all covariates except 'GenMetS'
        X_train = X_train.drop(columns=['GenMetS'])
        X_test = X_test.drop(columns=['GenMetS'])
    else:
        raise ValueError(f'Unknown input_config: {input_config}')
    # Use the xgb_best function to find the best hyperparameters
    xgb_best_model = xgb_best(X_train, y_train)
    # Create a bootstrap sample from the training set with equal number of cases and controls
    # Train the xgb model on the bootstrap sample using the best hyperparameters
    xgb_best_model.fit(X_train, y_train)
    
    # Evaluate the model on the testing set to compute the ROC curve metrics
    y_pred_prob = xgb_best_model.predict_proba(X_test)[:,1]
    fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
    roc_auc = auc(fpr, tpr)
    c_index_value = c_index(y_test.values, y_pred_prob)
    mean_prob_pos = y_pred_prob[y_test == 1].mean()
    mean_prob_neg = y_pred_prob[y_test == 0].mean()
    discrimination_slope = mean_prob_pos - mean_prob_neg
    
    # Calculate Calibration Slope
    calibration_model = LogisticRegression()
    calibration_model.fit(y_pred_prob.reshape(-1, 1), y_test)
    calibration_slope = calibration_model.coef_[0][0]
    # Compute SHAP values for feature importance
    explainer = shap.Explainer(xgb_best_model, X_train)
    shap_values = explainer.shap_values(X_test, check_additivity=False)
    covariates = list(X_test.columns)
    
    if input_config != 'GenMetS' and encoding_method == 'dummy':
        # Sum the absolute SHAP values for each level of Education
        # Education_2.0, Education_3.0, ..., Education_6.0 are the one-hot encoded columns for Education.
        education_shap_values = np.abs(shap_values[:, covariates.index('Education_2.0'): covariates.index('Education_6.0')+1]).sum(axis=1)
        
        # Replace one row of shap values for Education with the summed values
        shap_values[:, covariates.index('Education_2.0')] = education_shap_values
        
        # Delete the other rows of shap values for Education
        shap_values = np.delete(shap_values, np.s_[covariates.index('Education_3.0'): covariates.index('Education_6.0')+1], axis=1)
        
        # Adjust covariates list to match the modified shap_values
        covariates = [col for col in covariates if not col.startswith('Education_')] + ['Education']
    
    # Now compute mean absolute SHAP values
    vals = np.abs(shap_values).mean(0)
    
    # Create DataFrame for feature importance
    feature_importance = pd.DataFrame(list(zip(covariates, vals)), columns=['col_name', 'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'], ascending=False, inplace=True)
    
    return fpr, tpr, roc_auc, feature_importance, c_index_value, discrimination_slope, calibration_slope  # Add discrimination_slope and calibration_slope to the returned values

covariates_only = [cov for cov in covariates if cov != 'GenMetS']
# Initialize dictionaries to collect the results
model_specs = {
    'GenMetS + covariates': covariates,
    'GenMetS only': covariates_only,
    'GenMetS + age': covariates
}

model_cofig = {
    'GenMetS + covariates': 'cov + GenMetS',
    'GenMetS only': 'GenMetS',
    'GenMetS + age': 'GenMetS + Age'
}

auc_results = {
    'model': [],
    'bootstrap_iteration': [],
    'fpr': [],
    'tpr': [],
    'auc': []
}

feature_importance_results = {
    'model': [],
    'bootstrap_iteration': [],
    'feature': [],
    'feature_importance': []
}
c_index_results = {
    'model': [],
    'bootstrap_iteration': [],
    'c_index': []
}
discrimination_slope_results = {
    'model': [],
    'bootstrap_iteration': [],
    'discrimination_slope': []
}

calibration_slope_results = {
    'model': [],
    'bootstrap_iteration': [],
    'calibration_slope': []
}


# Define the model specifications
fpr_fixed = np.linspace(0, 1, 100)
encoding_method = 'normalize'  ## Can switch to 'dummy'
for model_name, model_covariates in model_specs.items():
    for i in tqdm(range(100), desc=f'Bootstrap iterations for {model_name}'):
        # Train and evaluate the model
        X_train, y_train, X_test, y_test = prepare_data(study_data_v2, target_col, covariates, seed=i)
        fpr, tpr, roc_auc, feature_importance, c_index_value, discrimination_slope, calibration_slope = train_evaluate_model(
            X_train, y_train, X_test, y_test, model_covariates, input_config=model_cofig[model_name], encoding_method= encoding_method
        )
        tpr_interpolated = interp1d(fpr, tpr, kind="linear", fill_value="extrapolate", bounds_error=False)(fpr_fixed)
        c_index_results['model'].append(model_name)
        c_index_results['bootstrap_iteration'].append(i)
        c_index_results['c_index'].append(c_index_value)
        discrimination_slope_results['model'].append(model_name)
        discrimination_slope_results['bootstrap_iteration'].append(i)
        discrimination_slope_results['discrimination_slope'].append(discrimination_slope)
        
        calibration_slope_results['model'].append(model_name)
        calibration_slope_results['bootstrap_iteration'].append(i)
        calibration_slope_results['calibration_slope'].append(calibration_slope)

        # Collect AUC results
        for f, t in zip(fpr_fixed, tpr_interpolated):
            auc_results['model'].append(model_name)
            auc_results['bootstrap_iteration'].append(i)
            auc_results['fpr'].append(f)
            auc_results['tpr'].append(t)
            auc_results['auc'].append(roc_auc)
        
        # Collect feature importance results
        for feature, importance in zip(feature_importance['col_name'], feature_importance['feature_importance_vals']):
            feature_importance_results['model'].append(model_name)
            feature_importance_results['bootstrap_iteration'].append(i)
            feature_importance_results['feature'].append(feature)
            feature_importance_results['feature_importance'].append(importance)

# Convert the results to DataFrames
auc_results_df = pd.DataFrame(auc_results)
feature_importance_results_df = pd.DataFrame(feature_importance_results)
c_index_results_df = pd.DataFrame(c_index_results)
discrimination_slope_results_df = pd.DataFrame(discrimination_slope_results)
calibration_slope_results_df = pd.DataFrame(calibration_slope_results)
directory = f'./XGboost_100_three_models/{encoding_method}/'
auc_results_df.to_csv(f'{directory}auc_results_df.csv')
feature_importance_results_df.to_csv(f'{directory}feature_importance_results_df.csv')
c_index_results_df.to_csv(f'{directory}c_index_results_df.csv')
discrimination_slope_results_df.to_csv(f'{directory}discrimination_slope_results_df.csv')
calibration_slope_results_df.to_csv(f'{directory}calibration_slope_results_df.csv')