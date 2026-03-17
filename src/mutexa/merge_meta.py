import pandas as pd
import argparse
import time
from datetime import date, datetime, timedelta
import csv
import re
from tqdm import tqdm
import warnings
import numpy as np
from typing import Optional



CONTINENTS = {
    'Africa': ['Algeria', 'Angola', 'Benin', 'Botswana', 'Burkina Faso', 'Burundi', 'Cabo Verde', 'Cameroon', 'Central African Republic', 'Chad', 'Comoros', 'Democratic Republic of the Congo', 'Djibouti', 'Egypt', 'Equatorial Guinea', 'Eritrea', 'Eswatini (formerly Swaziland)', 'Ethiopia', 'Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Ivory Coast', 'Kenya', 'Lesotho', 'Liberia', 'Libya', 'Madagascar', 'Malawi', 'Mali', 'Mauritania', 'Mauritius', 'Morocco', 'Mozambique', 'Namibia', 'Niger', 'Nigeria', 'Rwanda', 'Sao Tome and Principe', 'Senegal', 'Seychelles', 'Sierra Leone', 'Somalia', 'South Africa', 'South Sudan', 'Sudan', 'Tanzania', 'Togo', 'Tunisia', 'Uganda', 'Zambia', 'Zimbabwe'] ,
    'Asia': ['Afghanistan', 'Armenia', 'Azerbaijan', 'Bahrain', 'Bangladesh', 'Bhutan', 'Brunei', 'Cambodia', 'China', 'Cyprus', 'Georgia', 'India', 'Indonesia', 'Iran', 'Iraq', 'Israel', 'Japan', 'Jordan', 'Kazakhstan', 'Kuwait', 'Kyrgyzstan', 'Laos', 'Lebanon', 'Malaysia', 'Maldives', 'Mongolia', 'Myanmar (Burma)', 'Nepal', 'North Korea', 'Oman', 'Pakistan', 'Palestine', 'Philippines', 'Qatar', 'Saudi Arabia', 'Singapore', 'South Korea', 'Sri Lanka', 'Syria', 'Tajikistan', 'Thailand', 'Timor-Leste', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Uzbekistan', 'Vietnam', 'Yemen'],
    'Europe': ['Albania', 'Andorra', 'Austria', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Italy', 'Kosovo', 'Latvia', 'Liechtenstein', 'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Monaco', 'Montenegro', 'Netherlands', 'North Macedonia', 'Norway', 'Poland', 'Portugal', 'Romania', 'Russia', 'San Marino', 'Serbia', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom', 'Vatican City'],
    'North America': ['USA','Antigua and Barbuda', 'Bahamas', 'Barbados', 'Belize', 'Canada', 'Costa Rica', 'Cuba', 'Dominica', 'Dominican Republic', 'El Salvador', 'Grenada', 'Guatemala', 'Haiti', 'Honduras', 'Jamaica', 'Mexico', 'Nicaragua', 'Panama', 'Saint Kitts and Nevis', 'Saint Lucia', 'Saint Vincent and the Grenadines', 'Trinidad and Tobago', 'United States'],
    'Oceania': ['Australia', 'Fiji', 'Kiribati', 'Marshall Islands', 'Micronesia', 'Nauru', 'New Zealand', 'Palau', 'Papua New Guinea', 'Samoa', 'Solomon Islands', 'Tonga', 'Tuvalu', 'Vanuatu'],
    'South America': ['Argentina', 'Bolivia', 'Brazil', 'Chile', 'Colombia', 'Ecuador', 'Guyana', 'Paraguay', 'Peru', 'Suriname', 'Uruguay', 'Venezuela']
}

def get_continent(country_name):
    for continent, countries in CONTINENTS.items():
        if country_name in countries:
            return continent
    return None  # Return None if the country is not found in any continent


def data_preprocessing(df):
    # change to column name to the desired one
    # change the date format to desired one
    # if the host value is absent fill it with unknown value to avoid error later
    col_dict = dict()
    for name in df.columns:
        if re.search(r'accession|id', name, re.IGNORECASE):
            sampleid = name
            sampleid1 = 'sampleID'
            col_dict = {sampleid: sampleid1}
        elif re.search(r'^host$', name, re.IGNORECASE):
            host = name
            host1 = 'host'
            df.loc[df[host].isnull(), host] = 'unknown'
            col_dict.update({host: host1})
        elif re.search(r'country', name, re.IGNORECASE):
            country = name
            country1 = 'country'
            col_dict.update({country: country1})
        elif re.search(r'continent', name, re.IGNORECASE):
            country = name
            country1 = 'continent'
            col_dict.update({country: country1})
        elif re.search(r'state|province', name, re.IGNORECASE):
            country = name
            country1 = 'State'
            col_dict.update({country: country1})
        elif re.search(r'location', name, re.IGNORECASE):
            location = name
            location1 = 'location'
            col_dict.update({location: location1})
        elif re.search(r'collection|date', name, re.IGNORECASE):
            collectDate = name
            # -------------------------------
            print("Checking the format of the date column is correct ... ")
            # Drop rows where the date column has an empty value
            df.dropna(subset=[collectDate], inplace=True)
            # check the format of date column
            pattern_day = re.compile(r'\d{4}-\d{2}-XX')
            # Create a mask for rows that meet the condition
            mask = df[collectDate].str.match(pattern_day, na=False)
            # Check if any row meets the condition
            if mask.any():
                # Prompt the user
                choice = input(
                    f"Do you want to replace XX with 01 for date {df[collectDate].loc[mask].iloc[0]}? (y/n): ").lower()
                if choice == 'y':
                    # Replace 'XX' with '01' in the selected row(s)
                    df.loc[mask, collectDate] = df.loc[mask, collectDate].str.replace('XX', '01')
                    print(f"Replaced XX with 01 for date: {df[collectDate].loc[mask].iloc[0]}")
                elif choice == 'n':
                    # Remove the row(s) that meet the condition
                    print(f"Removing row(s) with date: {df[collectDate].loc[mask].iloc[0]}")
                    df = df[~mask]
                else:
                    print("Error: Invalid choice. Please enter 'yes' or 'no.")
                    exit(1)

            pattern_month_day = re.compile(r'\d{4}-XX-XX')
            mask = df[collectDate].str.match(pattern_month_day, na=False)
            rows_to_remove = df[mask]
            # Remove the rows and inform the user
            if not rows_to_remove.empty:
                print(f"Removing rows with dates: {', '.join(rows_to_remove[collectDate])}")
                df = df[~mask]

            df = df.loc[df[collectDate].str.len() >= 8]
            df.loc[:, collectDate] = pd.to_datetime(df[collectDate], dayfirst=True).dt.strftime('%Y-%m-%d')
            collectDate1 = 'collectDate'
            col_dict.update({collectDate: collectDate1})
        elif re.search(r'^age| age', name, re.IGNORECASE):
            age = name
            age1 = 'age'
            col_dict.update({age: age1})
        elif re.search(r'name', name, re.IGNORECASE):
            virus_name = name
            virus_name1 = 'virusName'
            col_dict.update({virus_name: virus_name1})
        elif re.search(r'lineage', name, re.IGNORECASE):
            lineage = name
            lineage1 = 'lineage'
            col_dict.update({lineage: lineage1})
        elif re.search(r'status', name, re.IGNORECASE):
            status = name
            status1 = 'status'
            col_dict.update({status: status1})
        elif re.search(r'gender', name, re.IGNORECASE):
            gender = name
            gender1 = 'gender'
            col_dict.update({gender: gender1})
        else:
            df.drop([name], axis=1, inplace=True)
    # rename the columns to the new names by using dictionary
    df.rename(columns=col_dict, inplace=True)
    # add column WHO at the end of our dataset
    df.insert(loc=len(df.columns), column='WHO', value=None)
    # df = df[~df['collectDate'].isna()]
    df.fillna('unknown', inplace=True)
    return df, col_dict

def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)







def run_merge(
    start: str,
    end: str,
    meta: str,
    threshold: float,
    mut: str,
    prefix: str,
    days: int,
    category: str,
):
        # Check if startDate is provided
    if start:
        start = start
    else:
        # If startDate is not provided, set it to 12 months before today
        start = date.today() - timedelta(days=365)

    if end:
        end = end
    else:
        end = date.today()

    try:
        R2_threshold = float(threshold)
    except ValueError:
        R2_threshold = 0.0

    metaFile = meta
    mutFile = mut

    # Ensure minDays is set correctly
    if days is None or days == '':
        minDays = 14
    else:
        try:
            minDays = int(days)
        except ValueError:
            minDays = 14

    # Ensure userColumn is set correctly
    if category is None or category == '':
        userColumn = "country"
    else:
        userColumn = category

    # Ensure prefix is set correctly
    if prefix is None or prefix == '':
        prefix = "output"
    else:
        prefix = prefix

    # startDate = "2019-01-01"
    # endDate = "2021-02-01"
    # R2_threshold = 0
    # metaFile = "ExampleData/metadata.csv"
    # mutFile = "ExampleData/mutation.csv"
    # prefix = "example"
    # minDays = 28
    # userColumn = 'lineage'

    mutlist_path = 'outputs/' + str(prefix) + '_mutlist.csv'
    hapdict_path = 'outputs/' + str(prefix) + '_hapdict.csv'

    # Threholds for mutations
    # minDays = 14
    minSamples = 10
    cdc = 'inputs/CDC_class.csv'


    continent_names = list(CONTINENTS.keys())

    print("Data Pre-Processing .... ")
    start_datetime = datetime.strptime(str(start), "%Y-%m-%d")
    end_datetime = datetime.strptime(str(end), "%Y-%m-%d")
    date_generated = [start_datetime + timedelta(days=x) for x in range(0, (end_datetime - start_datetime).days)]
    time = []
    times2days = dict()
    for date in date_generated:
        times = date.strftime("%Y-%m-%d")
        time.append(times)

    i = 1
    for times in time:
        times2days[times] = i
        i += 1

    # print(times2days['2020-02-11'])
    cdc_dict = {}
    with open(cdc, 'r', newline='') as f:
        reader = csv.reader(f, delimiter=',')
        for cdc, lin in reader:
            cdc_dict[cdc] = lin

    # read meta file into csv file
    df_metadata = pd.read_csv(metaFile)
    # make uniform format for column name and date and host
    df_metadata, dict_cols = data_preprocessing(df_metadata)
    print("Data Pre-Processing is done.")
    # print(df_metadata.columns)
    for index, row in tqdm(df_metadata.iterrows(), leave=False, desc="Parsing samples"):
        if 'lineage' in df_metadata.columns:
            # if we have a lineage in our data here we match this with CDC format and based on fill the WHO column
            target = next((s for s in cdc_dict.values() if str(row['lineage']) in s), None)
            if not (target is None):
                df_metadata.at[index, 'WHO'] = str(list(cdc_dict.keys())[list(cdc_dict.values()).index(target)])
        if start_datetime.date() <= pd.Timestamp(df_metadata.at[index, 'collectDate']).date() < end_datetime.date():
            # if the date of sample collection is in our target range and the host is human we add column day and
            # country to our dataframe.
            df_metadata.at[index, 'Days'] = times2days[df_metadata.at[index, 'collectDate']]
            if 'location' in df_metadata.columns:
                # Check if 'location' column has information
                if pd.notna(row['location']):
                    parts = row['location'].split('/')
                    # Check if the first part is a continent
                    part0 = parts[0].strip()
                    if part0 in continent_names:
                        df_metadata.at[index, 'country'] = parts[1].strip()
                        df_metadata.at[index, 'state'] = parts[2].strip() if len(parts) > 2 else None
                        df_metadata.at[index, 'city'] = parts[3].strip() if len(parts) > 3 else None
                        df_metadata.at[index, 'continent'] = part0
                    else:
                        df_metadata.at[index, 'country'] = part0
                        df_metadata.at[index, 'state'] = parts[2].strip() if len(parts) > 2 else None
                        df_metadata.at[index, 'city'] = parts[3].strip() if len(parts) > 3 else None
                        df_metadata.at[index, 'continent'] = get_continent(part0)

    # df_metadata = df_metadata[df_metadata['host'].str.lower() == 'human']
    if 'host' in df_metadata.columns:
        df_metadata = df_metadata[df_metadata['host'].str.lower().isin(['human', 'homo sapiens'])]

    if 'country' in df_metadata.columns:
        country_counts = df_metadata.groupby('country').size().reset_index(name='SamplesCount')
        filtered_countries = country_counts[country_counts['SamplesCount'] >= minSamples]['country']
        df_metadata = df_metadata.merge(filtered_countries, on='country')
    print("Parsing samples is done.")

    multlist = pd.read_csv(mutFile, sep='[:]', engine='python', header=None)
    mutations = list(multlist[0])

    hap_df = pd.read_csv(hapdict_path, header=None)
    hap_df.rename(columns={hap_df.columns[0]: 'sampleID'}, inplace=True)
    rows = hap_df[hap_df['sampleID'].str.contains('\.')]
    rows['sampleID'] = rows['sampleID'].str.split('\.').str[0]
    hap_df.update(rows)

    # change the column name based on mutaion name
    for i in range(1, len(hap_df.columns)):
        hap_df.rename(columns={hap_df.columns[i]: mutations[i - 1]}, inplace=True)


    # merge hap_df and metadata based on sampleID
    mutData = pd.merge(df_metadata, hap_df, on='sampleID')
    mutData.to_csv('outputs/' + prefix + "_ratioData.csv", index=False)

    # categorize data based on userColumn
    print("Categorize data based on the " + userColumn + '...')
    selectedMutations = {}
    if userColumn.lower() == 'continent':
        supergroup = mutData.continent.unique()
    if userColumn.lower() == 'country':
        supergroup = mutData.country.unique()
    if userColumn.lower() == 'lineage':
        supergroup = mutData.lineage.unique()
    ratio_out = 'outputs/' + prefix + '_' + userColumn + 'Ratios_' + str(R2_threshold) + '.csv'
    for item in tqdm(supergroup, 'R2 ratio analysis'):
        selectedMutations[item] = {}
        if userColumn.lower() == 'continent':
            selectedData = mutData[mutData.continent.eq(item)]
        if userColumn.lower() == 'country':
            selectedData = mutData[mutData.country.eq(item)]
        if userColumn.lower() == 'lineage':
            selectedData = mutData[mutData.lineage.eq(item)]
        for m in mutations:
            subset = selectedData[['Days', m]]
            subset.columns = ['Days', 'Mutation']
            SNP_present = subset[subset.Mutation.eq(1)]
            SNP_not_present = subset[subset.Mutation.eq(0)]
            SNPSamples = SNP_present.shape[0]
            RefSamples = SNP_not_present.shape[0]
            useOldMethod = True
            if SNPSamples > minSamples and RefSamples > minSamples and (not SNP_present['Days'].isna().all()):
                SNP_present = SNP_present.replace(np.nan, 0.00)
                SNP_not_present = SNP_not_present.replace(np.nan, 0.00)
                SNPDays = SNP_present["Days"].astype('int').max() - SNP_present['Days'].astype(
                    'int').min() + 1  # days between first and last SNP sample
                RefDays = SNP_not_present['Days'].astype('int').max() - SNP_not_present['Days'].astype(
                    'int').min() + 1  # days between first and last Ref sample
                if SNPDays >= minDays and RefDays >= minDays:
                    useOldMethod = False
                    SNPRate = SNPSamples / SNPDays
                    RefRate = RefSamples / RefDays
                    Ratio = SNPRate / RefRate
                    selectedMutations[item][m] = Ratio
            if useOldMethod:
                selectedMutations[item][m] = 0
    MutationRatios = pd.DataFrame.from_dict(selectedMutations, orient="index")

    output = MutationRatios[MutationRatios >= R2_threshold].dropna(thresh=1)
    output.fillna(0, inplace=True)
    # print(output)
    output.to_csv(ratio_out)
    print("Categorize data based on the " + userColumn + " is done.")

    print("Done!")


warnings.filterwarnings("ignore")


from typing import Optional


    


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Merge metadata and prepare MuTEXA tables")
    parser.add_argument("-s", "--start", required=True)
    parser.add_argument("-e", "--end", required=True)
    parser.add_argument("-m", "--meta", required=True)
    parser.add_argument("-t", "--thresh", default=0, type=float)
    parser.add_argument("-u", "--mut", required=True)
    parser.add_argument("-p", "--prefix", default="output")
    parser.add_argument("-d", "--days", default=14, type=int)
    parser.add_argument("-c", "--cat", required=True)

    args = parser.parse_args()

    run_merge(
        start=args.start,
        end=args.end,
        meta=args.meta,
        threshold=args.thresh,
        mut=args.mut,
        prefix=args.prefix,
        days=args.days,
        category=args.cat,
    )


if __name__ == "__main__":
    main()
