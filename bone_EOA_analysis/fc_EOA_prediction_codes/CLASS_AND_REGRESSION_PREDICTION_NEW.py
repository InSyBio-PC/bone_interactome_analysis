import numpy as np
import pandas as pd
import joblib
from scipy import stats
from sklearn.impute import KNNImputer
from sklearn.preprocessing import MinMaxScaler

# Load classification model
try:
    with open('/path/00001finalSingleModel.pkl.z', 'rb') as file:
        classification_model = joblib.load(file)
    print("Classification model loaded successfully.")
except Exception as e:
    print(f"Error loading the classification model: {e}")

# Load feature selection file
filter_file_path = '/path/clf_best_models/features_selected_1.csv'
filter_df = pd.read_csv(filter_file_path, index_col=0)
selected_features = filter_df.columns[filter_df.iloc[0] == 1].tolist()

# Load dataset
df = pd.read_csv("/path/ppis_with_features.csv")
print(df)
dfA = df.copy()

df = df.drop(columns=['uidA', 'uidB', 'protein_accession_A', 'protein_accession_B', 'seq_A', 'seq_B'])

# Handle NaN values
threshold = len(df) * 0.9
columns_to_replace = df.columns[df.isnull().sum() > threshold]
df[columns_to_replace] = df[columns_to_replace].fillna(0)

# Impute missing values and scale data
imputer = KNNImputer(n_neighbors=5)
df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns)
scaler = MinMaxScaler()
df_scaled = pd.DataFrame(scaler.fit_transform(df_imputed), columns=df.columns)
df = df_scaled

df1 = df[selected_features]

# Classification prediction
predicted_classes = classification_model.predict(df1)
predicted_probabilities = classification_model.predict_proba(df1)
predicted_class_probs = predicted_probabilities[np.arange(len(predicted_probabilities)), predicted_classes]

# Store classification results
result_df = pd.DataFrame({'Predicted Classes': predicted_classes,
                          'Probability Score': predicted_class_probs})
result_df = pd.concat([result_df, df], axis=1)
result_df['uidA'] = dfA['uidA']
result_df['uidB'] = dfA['uidB']
result_df.to_csv('/path/ppis_NEW_predictions.csv', index=False)

# Filter and save positive class predictions
filtered_df = result_df[result_df['Predicted Classes'] == 1]
filtered_df.to_csv('/path/ppi_combs_NEW_positives.csv', index=False)

# Regression Model Loading and Prediction
def load_model(model_path):
    with open(model_path, 'rb') as file:
        return joblib.load(file)

# Load regression models
model_info = [
    {'path': f'/path/Datasets/models/Output_dg_dataset_NEW_2/classification_models/000{str(i).zfill(2)}finalSingleModel.pkl.z', 'feature_row': i - 1}
    for i in range(1, 25) if i not in [12]
]
models = [load_model(info['path']) for info in model_info]

# Load feature selection data
features_1 = pd.read_csv('/path/Datasets/models/Output_dg_dataset_NEW_2/feature_selection/features_FinalFront1.csv')

# Perform regression predictions
ensemble_predictions = np.zeros(len(df))
for i, info in enumerate(model_info):
    selected_columns = features_1.iloc[info['feature_row']] == 1
    selected_feature_names = features_1.columns[selected_columns]
    filtered_test_data = df[selected_feature_names]
    predictions = models[i].predict(filtered_test_data)
    ensemble_predictions += predictions

ensemble_predictions /= len(models)
result_df['Regression Value'] = ensemble_predictions

# Save final results with regression values
result_df.to_csv('/path/ALL_regression_predictions.csv', index=False)
