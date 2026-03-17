import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
from datetime import date, datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import numpy as np
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
import pathlib
import matplotlib.ticker as mtick
import warnings
import json


def initialization(csvfile):
    countryGrouping = pd.read_csv(csvfile)
    countryGrouping = countryGrouping[['Sub-region Name', 'Country or Area']]
    country2region = {}
    for index, row in countryGrouping.iterrows():
        country2region[row[1]] = row[0]
    country2region['USA'] = 'Northern America'
    country2region['Taiwan'] = 'Eastern Asia'
    country2region['Hong Kong'] = 'Eastern Asia'
    country2region['United Kingdom'] = 'Northern Europe'
    country2region['South Korea'] = 'Eastern Asia'
    country2region['Vietnam'] = 'South-eastern Asia'
    country2region['Iran'] = 'Southern Asia'
    country2region['Czech Republic'] = 'Eastern Europe'
    country2region['Russia'] = 'Eastern Europe'
    country2region['Crimea'] = 'Eastern Europe'
    country2region['Brunei'] = 'South-eastern Asia'
    country2region['\u200eRomania'] = 'Eastern Europe'
    country2region['Venezuela'] = 'Latin America'
    country2region['Moldova'] = 'Eastern Europe'
    country2region['Reunion'] = 'Eastern Africa'
    country2region['Curacao'] = 'Latin America'
    country2region['Republic of Congo'] = 'Sub-Saharan Africa'
    country2region['Palestine'] = 'Western Asia'
    country2region['Saint Barthélemy'] = 'Latin America and the Caribbean'
    country2region['Saint Martin'] = 'Latin America and the Caribbean'
    country2region['CotedIvoire'] = 'Sub-Saharan Africa'
    country2region['Czech Repubic'] = 'Eastern Europe'
    country2region['St Eustatius'] = 'Latin America and the Caribbean'
    countryGrouping = pd.DataFrame.from_dict(country2region, orient='index')
    countryGrouping.columns = ['Region']
    countryGrouping.index.names = ['Country']
    return country2region, countryGrouping


def openRatioFile(prefix,userColumn):
    ratiodata = pd.read_csv('outputs/' + prefix + "_ratioData.csv")
    sortedData = ratiodata.sort_values([userColumn, "Days"], ascending=[True, True])
    sortedData["Date"] = pd.to_datetime(sortedData["collectDate"])
    date_series = pd.to_datetime(sortedData["Date"])
    date_index = pd.DatetimeIndex(date_series.values)
    sortedData.set_index("Date", inplace=True)
    sortedData = sortedData.set_index(date_index)
    sortedData = sortedData.rename_axis("Date")
    return sortedData


def openR2file(threshold, userColumn, prefix):
    R2_data = pd.read_csv('outputs/' + prefix + "_" + userColumn + "Ratios_" + str(threshold) + ".csv", dtype=object)
    R2_data.columns.values[0] = userColumn
    return R2_data


def add_strain(sortedData, mutations):
    omicron = r'B.1.1.529'
    BA_1 = r'BA.1|BA.1.*'
    BA_2 = r'BA.2|BA.2.*'
    BA_3 = r'BA.3|BA.3.*'
    BA_4 = r'BA.4|BA.4.*'
    BA_5 = r'BA.5|BA.5.*'
    delta = r'B.1.617.2|AY.*'
    alpha = r'^B\.1\.1\.7$|^Q.*'
    beta = r'B.1.351|B.1.351.*'
    gamma = r'^P.1|^P.1.*'
    mu = r'B.1.621|B.1.621.1'
    eta = 'B.1.525'
    iota = 'B.1.526'
    kappa = 'B.1.617.1'
    covlambda = 'C.37'
    c12 = 'C.1.2'
    epsilon = r'B.1.427|B.1.429'
    zeta = 'P.2'
    theta = 'P.3'
    xe = 'XE'
    xd = 'XD'

    sortedData[mutations] = sortedData[mutations].apply(pd.to_numeric)
    sortedData["lineage"] = sortedData["lineage"].fillna('None')

    sortedData['WHO'] = ''
    conditions = [
        (sortedData['lineage'].str.contains(alpha)),
        (sortedData['lineage'].str.contains(delta)),
        (sortedData['lineage'].str.contains(beta)),
        (sortedData['lineage'].str.contains(gamma)),
        (sortedData['lineage'].str.contains(omicron)),
        (sortedData['lineage'].str.contains(BA_1)),
        (sortedData['lineage'].str.contains(BA_2)),
        (sortedData['lineage'].str.contains(BA_3)),
        (sortedData['lineage'].str.contains(BA_4)),
        (sortedData['lineage'].str.contains(BA_5)),
        (sortedData['lineage'].str.contains(mu)),
        (sortedData['lineage'].str.contains(eta)),
        (sortedData['lineage'].str.contains(iota)),
        (sortedData['lineage'].str.contains(kappa)),
        (sortedData['lineage'].str.contains(covlambda)),
        (sortedData['lineage'].str.contains(c12)),
        (sortedData['lineage'].str.contains(epsilon)),
        (sortedData['lineage'].str.contains(zeta)),
        (sortedData['lineage'].str.contains(theta)),
        (sortedData['lineage'].str.contains(xe)),
        (sortedData['lineage'].str.contains(xd)),
        (sortedData['lineage'].str.contains('None'))
    ]
    choices = ['Alpha', 'Delta', 'Beta', 'Gamma', 'Omicron', 'BA.1', 'BA.2', 'BA.3', 'BA.4', 'BA.5', 'Mu', 'Eta',
               'Iota', 'Kappa', 'Lambda', 'c12', 'Epsilon', 'Zeta', 'Theta', 'XE', 'XD', 'None']
    sortedData['WHO'] = np.select(conditions, choices, default='Other')

    cols = sortedData.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    sortedData = sortedData[cols]
    return sortedData



def run_plots(
    start: str,
    end: str,
    threshold: float,
    mut: str,
    prefix: str,
    category: str,
    days: int,
) -> None:
    mutFile = mut
    threshold = float(threshold)
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

    # Ensure userColumn is set correctly
    if category is None or category == '':
        userColumn = "country"
    else:
        userColumn = category
    userColumn2 = userColumn.replace('"', '')

    # Ensure minDays is set correctly
    if days is None or days == '':
        Days = 14
    else:
        try:
            Days = int(days)
        except ValueError:
            Days = 14

    period = str(Days)+"D"

    # Ensure prefix is set correctly
    if prefix is None or prefix == '':
        prefix = "output"
    else:
        prefix = prefix

    # start = "2020-01-01"
    # end = "2021-10-04"
    # threshold = 0.0
    # mutFile = "ExampleData/mutation.csv"
    # prefix = ""
    # userColumn = "country"
    # userColumn2 = 'country'
    # Days = 14
    # period = str(Days)+"D"

    # print(countryGrouping)
    multlist = pd.read_csv(mutFile, sep='[:]', engine='python', header=None)
    mutations = list(multlist[0])
    ratioData = openRatioFile(prefix,userColumn)

    country2region, countryGrouping = initialization("inputs/UNSD — Methodology.csv")

    # ratioData = ratioData.rename(columns=name_map)
    cdata = openR2file(threshold, userColumn,prefix)
    allgroups = list(cdata[userColumn].unique())
    heatmapmatrix = {}

    for m in tqdm(mutations, desc='Calculating mutations rolling average'):
        print(m)
        ref = m
        subset = cdata[[userColumn, m]]
        subset.set_index(userColumn, inplace=True)
        subset.columns = ["Mutation"]
        allgroups = list(cdata[userColumn].unique())
        groupsOfInterest = list(subset[subset.Mutation.astype(float) > threshold].index)
        submut = ratioData[[userColumn, m, ref]]
        submut.columns = [userColumn, "Mutation", "Ref"]
        submut2 = submut.rename(columns={userColumn: userColumn, m: "Mutation", ref: "Ref"})
        countlist = ["Mutation", "Ref"]
        submut2[countlist] = submut2[countlist].astype(str).astype(float)
        submut2["Mutation"] = np.where(submut2["Mutation"].isin([1]), submut2["Mutation"], 0)
        submut2["Ref"] = np.where(submut2["Ref"].isin([0]), submut2["Ref"], np.nan)
        submut2["Ref"] = submut2["Ref"].map({np.nan: 0, 0: 1})
        countlist = ["Mutation", "Ref"]
        fill = [userColumn, submut2.index]
        # sum per day
        idx = pd.date_range(start, end, freq="1D")
        # Count per day mut/ref
        group = submut2.groupby(fill)[countlist].sum()
        group = group.reset_index()
        # print(group.shape)

        # fill in dates from Day 1 (2020-01-01) to end date
        group2 = (
            group.set_index("Date").groupby(userColumn).apply(lambda d: d.reindex(pd.date_range(start, end, freq='D')))
            .drop(userColumn, axis=1).reset_index(userColumn).fillna(np.nan))
        group2.rename_axis("Date", inplace=True)
        group2[countlist] = group2[countlist].astype(str).astype(float)
        # Rolling sum
        rolled = group2.groupby(userColumn)[countlist].rolling(period).sum()
        rolled["ratio"] = rolled["Mutation"] / rolled["Ref"]  # Calculate mutation ratios
        rolled2 = rolled[["ratio"]]
        rolled3 = rolled2.reset_index(userColumn)
        # mutdict = dict.fromkeys(mutations)
        allratios = {k: v for k, v in rolled3.groupby(userColumn)["ratio"]}
        mutdict = dict.fromkeys(mutations)
        # print(mutdict)
        cdict = {}
        # Check for infinity values in the 'ratio' column
        is_inf = np.isinf(rolled3['ratio'])
        # Countries with ratios over threshold
        # thres_pass = rolled3.loc[rolled3['ratio'] >= threshold]
        thres_pass = rolled3[~is_inf & (rolled3['ratio'] >= threshold)]
        cthres = list(thres_pass[userColumn2].unique())
        allgroupsofinterest = groupsOfInterest + list(set(cthres) - set(groupsOfInterest))

        # print(m,thres_pass['ratio'].max())
        # if thres_pass['ratio'].max() > 0: # xCarol -Save mutations to plot where max ratio is more than 0
        #     heatmap_to_plot_list.append(m) #Line 228 should work to do the same thing (tested on Flu and outputs same plots)

        for group in allgroupsofinterest:
            if np.nanmax(
                    rolled3.ratio) >= threshold:  # We should make it > to avoid plotting the zero ratios, but it gives me error[cannot do slice indexing on RangeIndex with these indexers [2022-01-01] of type str]

                cdict = dict.fromkeys(allgroupsofinterest)
                keylist = cdict.keys()
                for key in cdict:
                    cdict[key] = allratios[key]
            heatmapmatrix[m] = cdict
        mutdict.update(heatmapmatrix)

    # Plot heatmaps with countries above threshold
    mapp = mutdict
    for m in mapp:
        if not (mapp[m] is None):
            heatmaptoshow = pd.DataFrame.from_dict(mapp[m])
            # if not heatmaptoshow.empty:
            heatmaptoshow = heatmaptoshow[start:end]
            # if not heatmaptoshow.empty:
            results = heatmaptoshow.T
            results.reset_index()

            plotly_dict = {# March 2026
                    "type": "heatmap",
                    "x": results.columns.strftime('%Y-%m-%d').tolist(),        # dates
                    "y": results.index.tolist(),          # IIb, IIa, Ia etc
                    "z": results.values.tolist(),         # shape (len(y), len(x)) ✓
                    "showscale": True}  
            with open("outputs/"+prefix +"_" + m + str(threshold)+'_heatmapmatrix.json', 'w') as f:
                json.dump(plotly_dict, f)
            results = results.rename_axis(userColumn)
            
            if 'country' in results.columns:
                results = pd.merge(results, countryGrouping, how="inner", left_index=True, right_index=True)
                results = results.sort_values(['Region', userColumn2])
                results = results.drop(['Region'], axis=1)
            else:
                pass

            result = results
            result = result.rename_axis(userColumn)
            result = result.T
            result = result.rename_axis("Date")
            
            if not result.empty:
                plt.figure(figsize=(15, len(result.columns) / 10 + 5))
                mask = result.T.isnull()
                result_filled = result.fillna(-1)
                show = sns.heatmap(result.T, vmin=0,vmax=2, rasterized=True) # added vmin=2
                show.set_facecolor("black")
                ticklabels = [result.index[int(tick)].strftime('%Y-%m-%d') for tick in show.get_xticks()]
                show.set_xticklabels(ticklabels)
                plt.title(m, size=12)
                plt.yticks(rotation=0)
                show.set(ylabel=None)
                show.grid(False)
                plt.subplots_adjust(bottom=0.3)
                # plt.savefig("outputs/Heatmaps_mut_" + m + ".jpeg")
                plt.savefig("outputs/" + prefix + "_Heatmaps_mut:" + m + '_t:' + str(threshold) + ".jpeg")

            else:
                print(m)
                print('{}{}{}{}'.format("No",userColumn, "with ratio above ", threshold, " for this SNP\n"))
    # print("Saving heatmaps for: ", heatmap_to_plot_list)
    ratioData[mutations] = ratioData[mutations].apply(pd.to_numeric)
    # Save Strain counts for each mutation
    mutcounts_m = []
    for m in tqdm(mutations, desc='Saving mutation counts'):
        if len(ratioData[(ratioData[m] == 1) | (ratioData[m] == -1)]) > 0:
            desc = ratioData.reset_index()
            desc = ratioData[[userColumn2, m]] # Removed 'location'
            desc = desc.loc[desc[m] == 1]
            desc = desc.groupby([userColumn2]).describe()
            desc2 = desc[m]['count']
            desc2 = desc2.to_frame()
            desc2.columns.values[0] = m
            mutcounts_m.append(desc2)

    mutcounts_m2 = pd.concat(mutcounts_m, axis=1)
    mutcounts_m2 = mutcounts_m2.loc[:, ~mutcounts_m2.columns.duplicated()]
    mutcounts_m2 = mutcounts_m2.replace(np.nan, 0.00)
    #mutcounts_m2 = pd.DataFrame(mutcounts_m2, columns=[m]) # Carol - hashed
    topgroups = mutcounts_m2.sort_values(m, ascending=False).head(10).index.tolist()
    #mutcounts_m2.to_csv("outputs/"+ m+"mutcounts.csv") # Carol - for checking
    print(mutcounts_m2) # Show mutation counts for each group and mutation
    # Select mutations to plot rolling average if sum of mutations more than 10
    rolling_pass = mutcounts_m2.loc[:, (mutcounts_m2.sum() > 10)].columns.tolist() # Carol - Save mutation list for only those with sum >10

    counts = {}
    for m in tqdm(mutations):
        subset = ratioData[[userColumn, m]]
        subset.set_index(userColumn, inplace=True)
        subset.columns = ["Mutation"]

        submut = ratioData[[userColumn, m]]
        submut.columns = [userColumn, "Mutation"]
        submut2 = submut.rename(columns={userColumn: userColumn, m: "Mutation"})
        submut2["Mutation"] = np.where(submut2["Mutation"].isin([1]), submut2["Mutation"], 0)

        idx = pd.date_range(start, end, freq="1D")
        fill = [userColumn, submut2.index]
        group = submut2.groupby(fill)["Mutation"].sum()
        group = group.reset_index()
        group2 = (
            group.set_index("Date").groupby(userColumn).apply(lambda d: d.reindex(pd.date_range(start, end, freq='D')))
            .drop(userColumn, axis=1).reset_index(userColumn).fillna(np.nan))
        group2["Mutation"] = group2["Mutation"].astype(str).astype(float)
        rolled = group2.groupby(userColumn)["Mutation"].rolling(period).sum()
        rolled2 = rolled.reset_index(userColumn)
        rolled2 = rolled2.rename_axis("Date")

        countdict = dict.fromkeys(mutations)
        allcounts = {k: v for k, v in rolled2.groupby(userColumn)["Mutation"]}

        rolled3 = rolled2[[userColumn2, 'Mutation']]
        # appended_data.append(rolled3)
        rolled3.columns = [userColumn, m]
        # rolled3.to_csv("outputs/"+ m+"rolled.csv") # Carol - for checking
        for group in topgroups:
            cdict = dict.fromkeys(topgroups)
            keylist = cdict.keys()
            for key in cdict:
                cdict[key] = allcounts[key]
            counts[m] = cdict
        countdict.update(counts)
    ylabel = 'Rolling average ('+str(Days)+'-days)'

    print("Saving heatmaps for: ",len(heatmapmatrix.keys()),heatmapmatrix.keys())
    print("Saving rolling plots for: ", len(rolling_pass),rolling_pass)
    for m in rolling_pass: # Carol - added to use only mutation sums >10 to avoid plotting rolling averages that are too small
        linegraph = pd.DataFrame.from_dict(countdict[m])
        # linegraph.to_csv("outputs/"+ m+"linegraph.csv") # Carol - for checking
        # Get earlist date with mutation
        sub = ratioData[[userColumn2, m]]
        sub.columns = [userColumn2, "Mutation"]
        minD = sub.loc[sub["Mutation"].astype(float) > 0]
        mindate = minD.index.min().date()
        # Max rolling average
        max_value_column = pd.DataFrame(linegraph[1:].max())
        max_value = max_value_column.reset_index()
        max_value.columns = [userColumn2, 'value']
        maxgroup = max_value[max_value.value == max_value.value.max()]
        maxgroups = max_value.sort_values(by=['value'], ascending=False).head(10)
        maxgroups = maxgroups.loc[(maxgroups['value']!=0)]

        result = linegraph
        result = result.dropna(axis=1, how='all')  # Carol - Changed because it was just remove all columns?
        result = result[result.columns.intersection(list(maxgroups[userColumn2]))]  # Carol - only plot category in the top 10 that are not 0

        if not result.empty:
            sns.set_style("ticks")
            plt.figure(figsize=(12, 7))
            plt.xticks(rotation=60, horizontalalignment="right")
            plot = sns.lineplot(data=result,
                                dashes=False,
                                palette=sns.color_palette('tab10', n_colors=len(result.columns[:])), linewidth=2)
            plt.title(m, size=20)
            plot.set_ylabel(ylabel, fontsize=18)
            plot.set_xlabel("Date", fontsize=18)
            plot.tick_params(labelsize=16)
            date_form = DateFormatter("%Y-%m-%d")
            plot.xaxis.set_major_formatter(date_form)
            plot.set_xlim(datetime.datetime.strptime(start, "%Y-%m-%d"),datetime.datetime.strptime(end, "%Y-%m-%d"))# result.index.max().date()) # Carol - previously changed to same start, end for each mutation
            # plot.set_xlim(datetime.datetime.strptime(start, "%Y-%m-%d"), result.index.max().date())
            plot.xaxis.set_major_locator(mdates.WeekdayLocator(interval=5))
            plot.set_ylim(bottom=0)
            plot.axvspan(*mdates.datestr2num([str(mindate), str(mindate)]), color='red')
            plot.yaxis.set_major_locator(MaxNLocator(integer=True))
            plot.margins(y=0.05)
            maxdate = result[maxgroup.iloc[0][0]].idxmax()
            maxdate = maxdate.date()
            box = plot.get_position()
            plot.set_position([box.x0, box.y0 + box.height * 0.9,
                               box.width, box.height * 0.1])

            plt.tight_layout()
            plot.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
            plt.savefig("outputs/"+prefix+"_RollingAvg_mut:" + m + '_t:' + str(threshold) + ".jpeg", bbox_inches="tight", dpi=600)
            result.to_csv("outputs/" + prefix + "_RollingAvg_mut" +  m + "_t" + str(threshold) + ".csv") # March 2026
        else:
            print(m) # Carol - This line might be useless now after adding rolling_plots variable
            print("No observed mutations")



warnings.filterwarnings("ignore")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Create rolling average plots")
    parser.add_argument("-s", "--startDate")
    parser.add_argument("-e", "--endDate")
    parser.add_argument("-t", "--threshold", default=0)
    parser.add_argument("-u", "--mutFile", required=True)
    parser.add_argument("-p", "--prefix", default="output")
    parser.add_argument("-c", "--cat", default="country")
    parser.add_argument("-d", "--days", default=14)

    args = parser.parse_args()

    run_plots(
        start=args.startDate,
        end=args.endDate,
        threshold=args.threshold,
        mut=args.mutFile,
        prefix=args.prefix,
        category=args.cat,
        days=args.days,
    )

if __name__ == "__main__":
    main()