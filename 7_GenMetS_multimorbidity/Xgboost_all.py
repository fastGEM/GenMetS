# %%
from utils import *
import argparse

if __name__ == '__main__':
    ## Create the parser, add arguments and parse them for 'sex'
    parser = argparse.ArgumentParser(description="Parse user input for sex")
    # Add the 'sex' argument
    parser.add_argument(
        '--sex',
        type=str,
        choices=['Male', 'Female', 'Both'],
        help="Specify the sex of the individual (options: 'male', 'female', 'other')",
        default='Both'
    )

    # Parse the arguments
    args = parser.parse_args()
    
    target_col = 'multimorbidity'
    covariates = [
        "GenMetS", "TVhour", "Sleephour", "PhysicalActivity",
        "SocialEconomicalStatus", "Age", "Education"
    ]

    missing_covariates = ['GenMetS','Education']
    # Load the data
    df = pd.read_csv('./UKB_Asianwomenmen_s8792_mets_diseases_others.txt', sep='\t')
    ## Set the column Age as float
    df['Age'] = df['Age'].astype('float64')
    ## set eduction = NA when -3, -7
    df.loc[df.Education == -3, 'Education'] = np.nan
    df.loc[df.Education == -7, 'Education'] = np.nan
    female_df = df[df.sex == 'Female']
    female_df.reset_index(drop=True, inplace=True)
    male_df = df[df.sex == 'Male']
    male_df.reset_index(drop=True, inplace=True)
    if args.sex == 'Both':
        X_train, y_train, X_test, y_test = prepare_data(df, target_col, covariates, seed=42)
        print("Data prepared for Both sexes")
    else:
        per_sex_df = df[df.sex == args.sex]
        X_train, y_train, X_test, y_test = prepare_data(per_sex_df, target_col, covariates, seed=42)
        print(f"Data prepared for {args.sex}")
    if not os.path.exists(f'./{args.sex}/XGboost_100_three_models'):
        os.makedirs(f'./{args.sex}/XGboost_100_three_models')
    for encoding_method in ['dummy', 'normalize']:
        if not os.path.exists(f'./{args.sex}/XGboost_100_three_models/{encoding_method}'):
            os.makedirs(f'./{args.sex}/XGboost_100_three_models/{encoding_method}')

    model_cofig = {
        'GenMetS + covariates': 'cov + GenMetS',
        'GenMetS only': 'GenMetS',
        'GenMetS + age': 'GenMetS + Age',
        'Age Only': 'Age',
        'covariates only': 'cov'
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
    for model_name, model_covariates in model_cofig.items():
        for i in tqdm(range(100), desc=f'Bootstrap iterations for {model_name}'):
            # Train and evaluate the model
            fpr, tpr, roc_auc, feature_importance, c_index_value, discrimination_slope, calibration_slope = train_evaluate_model(
                X_train, y_train, X_test, y_test, input_config=model_cofig[model_name], encoding_method= encoding_method
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
    directory = f'./{args.sex}/XGboost_100_three_models/{encoding_method}/'
    auc_results_df.to_csv(f'{directory}auc_results_df.csv')
    feature_importance_results_df.to_csv(f'{directory}feature_importance_results_df.csv')
    c_index_results_df.to_csv(f'{directory}c_index_results_df.csv')
    discrimination_slope_results_df.to_csv(f'{directory}discrimination_slope_results_df.csv')
    calibration_slope_results_df.to_csv(f'{directory}calibration_slope_results_df.csv')

