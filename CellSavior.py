#!/usr/bin/env python3
#
# CellSavior: A tool for cell calling from target-based single-cell DNA sequencing data
# Author: Taehyeon Kim
# Date: 2026-04-10
# Version: 1.0.3
#
#--------------------------------#
# Import libraries
#--------------------------------#
import argparse
import pandas
import numpy as np
import umap
import hdbscan
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import pairwise_distances
from hdbscan.validity import validity_index
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from scipy.stats import rankdata
import scipy.signal
import pysam
#--------------------------------#


class CellSavior:
    def __init__(
        self, 
        input_barcode_distribution, 
        input_bam, 
        output_barcode_distribution, 
        output_bam, 
        output_prefix_for_plots, 
        min_reads_per_barcode=1.25, 
        min_covered_amplicons_per_barcode=0.40, 
        filter_amplicons=False,
        amplicon_variance_threshold=0.95,
        log_scaling=True,
        n_jobs=1, 
        thresholding_method="best", 
        completeness_threshold=0.8, 
        min_cell_barcodes_count=1000
    ):
        self.input_barcode_distribution = input_barcode_distribution
        self.input_bam = input_bam
        self.output_barcode_distribution = output_barcode_distribution
        self.output_bam = output_bam
        self.output_prefix_for_plots = output_prefix_for_plots
        self.min_reads_per_barcode = min_reads_per_barcode
        self.min_covered_amplicons_per_barcode = min_covered_amplicons_per_barcode
        self.amplicon_variance_threshold = amplicon_variance_threshold
        self.filter_amplicons = filter_amplicons
        self.log_scaling = log_scaling
        self.n_jobs = n_jobs
        self.thresholding_method = thresholding_method
        self.completeness_threshold = completeness_threshold
        self.min_cell_barcodes_count = min_cell_barcodes_count

    def select_cell_candidates(self, barcode_distribution):
        #STEP2
        N_amplicons = barcode_distribution.shape[1]
        total_reads = barcode_distribution.sum(axis=1)
        condition_1 = total_reads > self.min_reads_per_barcode * N_amplicons
        covered_amplicons = (barcode_distribution >= 1).sum(axis=1)
        condition_2 = covered_amplicons >= self.min_covered_amplicons_per_barcode * N_amplicons
        cell_candidates = barcode_distribution[condition_1 & condition_2]
        return cell_candidates

    def filter_amplicons(self, cell_candidates):
        #STEP2.1
        variance_per_amplicon = cell_candidates.var(axis=0)
        amplicon_filtered_cell_candidates = cell_candidates.loc[:, variance_per_amplicon < variance_per_amplicon.quantile(self.amplicon_variance_threshold)]
        return amplicon_filtered_cell_candidates

    def normalize_cell_candidates(self, cell_candidates):
        #STEP3
        if self.log_scaling == True:
            log_candidates = np.log1p(cell_candidates)
            normalized_candidates = log_candidates.sub(log_candidates.mean(axis=1), axis=0)
        else:
            normalized_candidates = cell_candidates.div(cell_candidates.sum(axis=1), axis=0)
        return normalized_candidates

    def UMAP(self, normalized_cell_candidates):
        #STEP4
        umap_metrics = [
            "euclidean", "manhattan", "chebyshev", "minkowski", "canberra", 
            "braycurtis", "seuclidean", "cosine", 
            "precomputed"
        ]
        precomputed_metrics = [
            "correlation", "sqeuclidean", "spearman"
        ]

        umap_reducers = {}
        for metric_to_use in umap_metrics:
            umap_reducers[metric_to_use] = umap.UMAP(
                densmap = False,
                n_neighbors = 20, 
                n_components = 2, 
                min_dist = 0, 
                spread = 1, 
                local_connectivity = 1, 
                metric = metric_to_use,
                n_jobs=self.n_jobs
            )
        
        embedded_umaps = {}
        for metric_to_use in umap_metrics:
            if metric_to_use == "precomputed":
                for metric in precomputed_metrics:
                    print(f"        Processing metric: {metric}")
                    if metric == "correlation":
                        distance_matrix = pairwise_distances(normalized_cell_candidates, metric="correlation", n_jobs=self.n_jobs)
                        correlation_matrix = 1 - distance_matrix.copy()
                    if metric == "sqeuclidean":
                        distance_matrix = pairwise_distances(normalized_cell_candidates, metric="sqeuclidean", n_jobs=self.n_jobs)
                    if metric == "spearman":
                        rank_cell_candidates = scipy.stats.rankdata(normalized_cell_candidates, axis=1)
                        distance_matrix = pairwise_distances(rank_cell_candidates, metric="correlation", n_jobs=self.n_jobs)

                    try:
                        embedded_umaps[metric] = umap_reducers[metric_to_use].fit_transform(distance_matrix)
                        print(f"        Metric {metric} embedded successfully.\n")
                    except Exception as e:
                        print(f"        WARNING::: embedding with metric {metric} failed: {e}\n")
                        continue
            else:
                try:
                    print(f"        Processing metric: {metric_to_use}")
                    embedded_umaps[metric_to_use] = umap_reducers[metric_to_use].fit_transform(normalized_cell_candidates)
                    print(f"        Metric {metric_to_use} embedded successfully.\n")
                except Exception as e:
                    print(f"        WARNING::: embedding with metric {metric_to_use} failed: {e}\n")
                    continue

        return embedded_umaps, correlation_matrix

    def HDBSCAN(self, embedded_umaps):
        #STEP5
        hdbscan_clusterer = hdbscan.HDBSCAN(
            metric = 'euclidean', 
            min_cluster_size = 15, 
            min_samples = 15, 
            cluster_selection_epsilon = 0.01, 
            cluster_selection_method = 'eom', 
            alpha = 1.0, 
            core_dist_n_jobs=self.n_jobs
        )

        dbcv_score = {}
        clustered_df = {}
        succeded_metrics = list(embedded_umaps.keys())

        for metric in succeded_metrics:
            print(f"        Processing metric: {metric}")
            hdbscan_clusterer.fit(embedded_umaps[metric])
            labels = hdbscan_clusterer.fit_predict(embedded_umaps[metric])

            dbcv_score[metric] = validity_index(
                np.array(embedded_umaps[metric], dtype = np.float64), 
                labels
            )
            clustered_df[metric] = pandas.DataFrame({
                "UMAP_1": embedded_umaps[metric][:, 0], 
                "UMAP_2": embedded_umaps[metric][:, 1], 
                "cluster_id": hdbscan_clusterer.labels_
            })
            print(f"        Metric {metric} clustered successfully.\n")
        
        print("        Calculating validity scores...")
        valid_metrics = {metric: np.exp(score) for metric, score in dbcv_score.items()}
        valid_metrics = {k: v/sum(valid_metrics.values()) for k, v in valid_metrics.items()}
        valid_metrics = dict(sorted(valid_metrics.items(), key=lambda x: x[1], reverse=True))
        metric_values_df = pandas.DataFrame({
            "metric": succeded_metrics,
            "dbcv_score": [dbcv_score[m] for m in succeded_metrics],
            "weight": [valid_metrics[m] for m in succeded_metrics]
        })
        print("        Validity scores calculated successfully.\n")

        return clustered_df, metric_values_df

    def calculate_validity(self, cell_candidates, normalized_cell_candidates, correlation_matrix, clustered_df, metric_values_df):
        #STEP6
        N_amplicons = len(normalized_cell_candidates.columns)

        correlation_coverages_df = pandas.DataFrame({
            "total_depth": np.log10(cell_candidates.sum(axis = 1) / N_amplicons), 
            "R-squared": (correlation_matrix**2).mean(axis = 1), 
            'ado': (cell_candidates == 0).sum(axis = 1) / N_amplicons
        })

        for metric in metric_values_df["metric"].values.tolist():
            correlation_coverages_df[metric] = clustered_df[metric]['cluster_id'].values.tolist()

        validation_df = correlation_coverages_df.copy()
        cluster_validities_by_metric = {}
        for metric in metric_values_df["metric"].values.tolist():
            print(f"        Processing metric: {metric}")
            print(f"            1.Feature extraction of clusters...")
            clusters = np.array(sorted(validation_df[metric].unique()))
            clusters = clusters[clusters != -1]

            cluster_values = {"total_depth":[], "R-squared":[], "ado":[], "density":[]}
            for cluster_id in clusters:
                cluster_info = validation_df[validation_df[metric] == cluster_id].iloc[:, 0:3]

                cluster_values['total_depth'].append(np.quantile(cluster_info['total_depth'], 0.10))
                cluster_values['R-squared'].append(np.quantile(cluster_info['R-squared'], 0.90))
                cluster_values['ado'].append(np.quantile(cluster_info['ado'], 0.90)**1.5)

                k=min(15, len(cluster_info)-1)
                knn = NearestNeighbors(n_neighbors=k+1).fit(cluster_info.iloc[:, 0:2])
                distances, _ = knn.kneighbors(cluster_info.iloc[:, 0:2])
                cluster_values['density'].append( 1.0 / np.quantile(distances[:, 1:].mean(axis = 1), 0.10) )

            cluster_validity = pandas.DataFrame(cluster_values)
            density_rank = rankdata(cluster_validity['density'], method='average')
            cluster_validity['density'] = density_rank / density_rank.max()
            print(f"            Feature extraction of clusters completed successfully.\n")

            print(f"            2.Binarizing clusters...")
            fmin, fmax = cluster_validity.min(axis = 0), cluster_validity.max(axis = 0)
            init_centroids = np.vstack([
                np.array([fmin.iloc[0], fmax.iloc[1], fmax.iloc[2], fmax.iloc[3]]),
                np.array([fmax.iloc[0], fmin.iloc[1], fmin.iloc[2], fmin.iloc[3]])
            ])

            kmeans = KMeans(n_clusters=2, init=init_centroids, n_init=1)
            cluster_validity['kmeans'] = kmeans.fit_predict(cluster_validity)

            cluster_1 = cluster_validity[cluster_validity['kmeans'] == 0]
            mean_values_1 = cluster_1.mean(axis = 0)

            cluster_2 = cluster_validity[cluster_validity['kmeans'] == 1]
            mean_values_2 = cluster_2.mean(axis = 0)
            print(f"            Binarizing clusters completed successfully.\n")


            print(f"            3.Checking validity of binarized groups...")
            checkpoints = {
                'total_depth': mean_values_1['total_depth'] < mean_values_2['total_depth'], 
                'R-squared': mean_values_1['R-squared'] > mean_values_2['R-squared'], 
                'ado': mean_values_1['ado'] > mean_values_2['ado'], 
                'density': mean_values_1['density'] > mean_values_2['density']
            }

            cluster_validity['validity'] = ""
            if sum(checkpoints.values()) >= 3:
                cluster_validity.loc[cluster_validity['kmeans'] == 0, 'validity'] = "invalid"
                cluster_validity.loc[cluster_validity['kmeans'] == 1, 'validity'] = "valid"
            elif sum(checkpoints.values()) <= 1:
                cluster_validity.loc[cluster_validity['kmeans'] == 0, 'validity'] = "valid"
                cluster_validity.loc[cluster_validity['kmeans'] == 1, 'validity'] = "invalid"
            else:
                cluster_validity.loc[cluster_validity['kmeans'] == 0, 'validity'] = "ambiguous"
                cluster_validity.loc[cluster_validity['kmeans'] == 1, 'validity'] = "ambiguous"

            cluster_validities_by_metric[metric] = cluster_validity
            print(f"            Validity of binarized groups checked successfully.\n")

        cluster_map = {
            metric: {
                cluster_id: 
                    1 if validity == "valid" 
                    else 0.5 if validity == "ambiguous" 
                    else 0 
                for cluster_id, validity in cluster_validity['validity'].items()
            }
            for metric, cluster_validity in cluster_validities_by_metric.items()
        }

        mapped_df = validation_df.iloc[:, 0:2].copy()
        for metric, mapping in cluster_map.items():
            mapped_df[metric] = validation_df[metric].map(mapping)

        mapped_df = mapped_df.fillna(0)
        mapped_df['total_score'] = sum(
            weight * mapped_df[metric]
            for metric, weight in zip(self.metric_values_df["metric"].values.tolist(), self.metric_values_df["weight"].values.tolist())
        )

        return cluster_map, mapped_df

    def validity_thresholding(self, mapped_df):
        #STEP7
        print("        STEP-7.1: Calculating KDE threshold...")
        kdeplot = sns.kdeplot(mapped_df['total_score']).lines[0]
        x_vals = kdeplot.get_xdata()
        y_vals = kdeplot.get_ydata()
        local_minima_indices = scipy.signal.argrelextrema(y_vals, np.less)[0]
        kde_threshold_df = pandas.DataFrame({
            "threshold_candidates": [x_vals[val] for val in local_minima_indices],
            "filtered_barcodes_count": [len(mapped_df[mapped_df['total_score'] > x_vals[val]]) for val in local_minima_indices]
        })
        print("        KDE threshold calculated successfully.\n")

        print("        STEP-7.2: Calculating elbow threshold...")
        th_ = [i/10000 for i in range(0, 10001)]
        counts_filtered = [len(mapped_df[mapped_df['total_score'] > th_[i]]) - len(mapped_df[mapped_df['total_score'] > th_[i+1]]) for i in range(0, 10000)]
        counts_filtered.insert(0, len(mapped_df) - len(mapped_df['total_score'] > 0))
        ellbow_threshold = (np.argmax(counts_filtered[1:5000]) + 1) / 10000
        print("        Elbow threshold calculated successfully.\n")

        print("        STEP-7.3: Completeness algorithm...")
        N_amplicons = self.barcode_distribution.shape[1]
        total_reads = self.barcode_distribution.sum(axis = 1)
        read_cutoff = N_amplicons * 8

        cell_candidates = self.barcode_distribution.loc[total_reads >= read_cutoff]
        panel_performance = 0.2 * cell_candidates.values.mean()
        amplicon_mean = cell_candidates.mean(axis = 0)

        good_amplicons = amplicon_mean[amplicon_mean >= panel_performance].index
        good_data = self.barcode_distribution[good_amplicons]

        completeness = (good_data > 0).sum(axis = 1) / len(good_amplicons)
        completeness_based_cell_barcodes = completeness[completeness >= self.completeness_threshold].index
        print("        Completeness algorithm completed successfully.\n")

        
        if self.thresholding_method == "best":
            print("        STEP-7.4: Searching for the best threshold...")
            kde_filtered_thresholds = kde_threshold_df.loc[
                    (kde_threshold_df['filtered_barcodes_count'] > self.min_cell_barcodes_count)
                    & (kde_threshold_df['filtered_barcodes_count'] < len(mapped_df)), 
                    'threshold_candidates'
                ].values.tolist()
            
            if len(kde_filtered_thresholds) > 0:
                final_threshold = np.min(kde_filtered_thresholds)
                cell_barcodes = mapped_df[mapped_df['total_score'] > final_threshold].index
                print("        The best threshold found successfully.\n")
                return final_threshold, cell_barcodes

            else:
                ellbow_based_cell_barcodes_count = len(mapped_df[mapped_df['total_score'] > ellbow_threshold].index)
                
                if ellbow_based_cell_barcodes_count > self.min_cell_barcodes_count:
                    final_threshold = ellbow_threshold
                    cell_barcodes = mapped_df[mapped_df['total_score'] > final_threshold].index
                    print("        The best threshold found successfully.\n")
                    return final_threshold, cell_barcodes

                else:
                    final_threshold = 0
                    cell_barcodes = completeness_based_cell_barcodes
                    print("        The best threshold found successfully.\n")
                    return final_threshold, cell_barcodes

        else:
            print("        STEP-7.4: Final thresholding...")
            if self.thresholding_method == "kde":
                final_threshold = np.min(kde_filtered_thresholds)
                cell_barcodes = mapped_df[mapped_df['total_score'] > final_threshold].index
                print("        Final thresholding completed successfully.\n")
                return final_threshold, cell_barcodes
            elif self.thresholding_method == "ellbow":
                final_threshold = ellbow_threshold
                cell_barcodes = mapped_df[mapped_df['total_score'] > final_threshold].index
                print("        Final thresholding completed successfully.\n")
                return final_threshold, cell_barcodes
            else:
                final_threshold = 0
                cell_barcodes = completeness_based_cell_barcodes
                print("        Final thresholding completed successfully.\n")
                return final_threshold, cell_barcodes

    def generate_barcode_distribution(self, cell_barcodes):
        output_data = self.barcode_distribution.loc[cell_barcodes]
        output_data.index.name = "cell_barcode"
        output_data.to_csv(self.output_barcode_distribution, sep="\t", header = True, index = True)

    def generate_cells_bam(self, cell_barcodes):
        input_bam = pysam.AlignmentFile(self.input_bam, "rb")
        header_dict = input_bam.header.to_dict()
        cell_barcodes = set(cell_barcodes)
        header_dict = input_bam.header.to_dict()
        if "RG" in header_dict:
            header_dict["RG"] = [rg for rg in header_dict["RG"] if rg["ID"] in cell_barcodes]
        output_bam = pysam.AlignmentFile(self.output_bam, "wb", header=header_dict)
        filtered_count = 0
        total_count = 0
        for read in input_bam.fetch(until_eof=True):
            total_count += 1
            if read.has_tag("RG"):
                barcode = read.get_tag("RG")
                if barcode in cell_barcodes:
                    output_bam.write(read)
                    filtered_count += 1
        input_bam.close()
        output_bam.close()

    def plot_clustering_results(self, clustered_data, metric_values_df):
        best_metric = metric_values_df.loc[metric_values_df['weight'].idxmax(), 'metric']
        grid = (-(-len(clustered_data.keys()) // 4), 4)
        fig, axes = plt.subplots(grid[0], grid[1], figsize=(20, 15), facecolor = "white", dpi = 100)
        k=0
        for i in range(grid[0]):
            for j in range(grid[1]):
                if k == len(metric_values_df['metric'].values.tolist()): break
                metric_to_use = metric_values_df['metric'].values.tolist()[k]
                plot_df = clustered_data[metric_to_use]
                sns.scatterplot(
                    data = plot_df[plot_df['cluster_id'] == -1],
                    x = "UMAP_1", 
                    y = "UMAP_2", 
                    color = "lightgray",
                    alpha = 0.3, 
                    s = 10, 
                    ax = axes[i, j], 
                    legend = False
                )
                sns.scatterplot(
                    data = plot_df[plot_df['cluster_id'] != -1],
                    x = "UMAP_1", 
                    y = "UMAP_2", 
                    hue = "cluster_id",
                    palette = "viridis",
                    alpha = 0.3, 
                    s = 10, 
                    ax = axes[i, j], 
                    legend = False
                )
                if metric_to_use == best_metric:
                    axes[i, j].set_title(f"{metric_to_use} :: {metric_values_df.loc[metric_values_df['metric'] == metric_to_use, 'weight'].values[0]:.2f}", color='blue')
                    for spine in axes[i, j].spines.values():
                        spine.set_edgecolor('blue')
                elif metric_to_use in metric_values_df['metric'].values.tolist():
                    axes[i, j].set_title(f"{metric_to_use} :: {metric_values_df.loc[metric_values_df['metric'] == metric_to_use, 'weight'].values[0]:.2f}", color='skyblue')
                    for spine in axes[i, j].spines.values():
                        spine.set_edgecolor('skyblue')
                else:
                    axes[i, j].set_title(f"{metric_to_use} :: {metric_values_df.loc[metric_values_df['metric'] == metric_to_use, 'weight'].values[0]:.2f}", color='lightgray')
                    for spine in axes[i, j].spines.values():
                        spine.set_edgecolor('lightgray')
                k+=1
        total_plots = len(metric_values_df['metric'].values.tolist())
        total_axes = grid[0] * grid[1]
        for idx in range(total_plots, total_axes):
            i = idx // 4
            j = idx % 4
            axes[i, j].axis("off")
        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        plt.savefig(f"{self.output_prefix_for_plots}.umap_hdbscan.results.png", facecolor = 'white', dpi = 350)

    def plot_final_correlation_coverage_plot(self, mapped_df, threshold = -1):
        plt.figure(dpi = 100, figsize = (10, 6), facecolor = "white")
        sns.scatterplot(
            data = mapped_df, 
            x = "total_depth", y = "R-squared", 
            hue = "total_score", 
            palette = "Blues"
        )
        plt.savefig(f"{self.output_prefix_for_plots}.correlation_coverage_plot_with_validity_score.png", facecolor = 'white', dpi = 350)
        plt.figure(dpi = 100, figsize = (10, 6), facecolor = "white")
        if threshold != -1:
            sns.scatterplot(
                data = mapped_df[mapped_df['total_score'] <= threshold], 
                x = "total_depth", y = "R-squared", 
                color = "lightgray", label = "invalid_barcodes"
            )
            sns.scatterplot(
                data = mapped_df[mapped_df['total_score'] > threshold], 
                x = "total_depth", y = "R-squared", 
                palette = "darkblue", label = "valid_cells"
            )
            plt.savefig(f"{self.output_prefix_for_plots}.correlation_coverage_plot_with_valid_cells.png", facecolor = 'white', dpi = 350)

    def plot_determination_of_validity_threshold(self, mapped_df, threshold):
        plt.figure(dpi = 100, figsize = (10, 6), facecolor = "white")
        sns.lineplot(
            x = [threshold / 10000 for threshold in range(0, 10001)], 
            y = [len(mapped_df[mapped_df['total_score'] > i / 10000]) for i in range(0, 10001)], 
            color = 'black'
        )
        plt.axvline(
            x = threshold, 
            color = 'darkred', 
            linestyle = '--'
        )
        plt.savefig(f"{self.output_prefix_for_plots}.filtered_count.lineplot.png", facecolor = 'white', dpi = 350)
        #
        plt.figure(dpi = 100, figsize = (10, 6), facecolor = "white")
        sns.kdeplot(mapped_df['total_score'])
        plt.axvline(
            x = threshold, 
            color = 'darkred', 
            linestyle = '--'
        )
        plt.savefig(f"{self.output_prefix_for_plots}.validity_score.kdeplot.png", facecolor = 'white', dpi = 350)



    def execute(self):
        print("    STEP-1: Loading barcode distribution matrix...")
        self.barcode_distribution = pandas.read_csv(
            self.input_barcode_distribution, 
            sep="\t"
            ).set_index("barcode")
        print("   Barcode distribution matrix loaded successfully.")
        print("\n\n")


        print("    STEP-2: Selecting cell candidates...")
        self.cell_candidates = self.select_cell_candidates(self.barcode_distribution)
        if self.filter_amplicons == True:
            print("        STEP-2.1: Filtering amplicons...")
            self.cell_candidates = self.filter_amplicons(self.cell_candidates)
            print("        Amplicons filtered successfully.")
        print("    Cell candidates selected successfully.")
        print("\n\n")


        print("    STEP-3: Normalizing cell candidates...")
        self.normalized_cell_candidates = self.normalize_cell_candidates(self.cell_candidates)
        print("    Cell candidates normalized successfully.")
        print("\n\n")


        print("    STEP-4: UMAP embedding...")
        self.embedded_umaps, self.correlation_matrix = self.UMAP(self.normalized_cell_candidates)
        print("    UMAP embedding completed successfully.")
        print("\n\n\n\n\n")


        print("    STEP-5: HDBSCAN clustering...")
        self.hdbscan_clusters, self.metric_values_df = self.HDBSCAN(self.embedded_umaps)
        print("    HDBSCAN clustering completed successfully.")
        print("\n\n")

        print("    STEP-6: Calculating validity scores...")
        self.cluster_map, self.mapped_df = self.calculate_validity(self.cell_candidates, self.normalized_cell_candidates, self.correlation_matrix, self.hdbscan_clusters, self.metric_values_df)
        print("    Validity scores calculated successfully.")
        print("\n\n")

        print("    STEP-7: Cell calling...")
        self.score_threshold, self.cell_barcodes = self.validity_thresholding(self.mapped_df)
        print("    Cell calling completed successfully.")
        print("\n\n")

        print("    STEP-8: Plotting results...")
        self.plot_clustering_results(self.hdbscan_clusters, self.metric_values_df)
        self.plot_final_correlation_coverage_plot(self.mapped_df, self.score_threshold)
        self.plot_determination_of_validity_threshold(self.mapped_df, self.score_threshold)
        print("    Results plotted successfully.")
        print("\n\n")

        print("    STEP-9: Generating barcode distribution...")
        self.generate_barcode_distribution(self.cell_barcodes)
        print("    Barcode distribution generated successfully.")
        print("\n\n")

        print("    STEP-10: Generating BAM file...")
        self.generate_cells_bam(self.cell_barcodes)
        print("    BAM file generated successfully.")
        print("\n\n")

#--------------------------------#
# Main method
#--------------------------------#
def parse_args():
    parser = argparse.ArgumentParser(
        description="CellSavior: A tool for cell calling from target-based single-cell DNA sequencing data"
        )
    
    parser.add_argument(
        "--input_barcode_distribution", 
        type=str, 
        required=True, 
        help="Input barcode distribution file path (required)"
        )
    parser.add_argument(
        "--input_bam", 
        type=str, 
        required=True, 
        help="Input BAM file path (required)"
        )
    parser.add_argument(
        "--output_barcode_distribution", 
        type=str, 
        required=True, 
        help="Output barcode distribution file path (required)"
        )
    parser.add_argument(
        "--output_bam", 
        type=str, 
        required=True, 
        help="Output BAM file path (required)"
        )
    parser.add_argument(
        "--output_prefix_for_plots", 
        type=str, 
        required=True, 
        help="Output prefix for plots (required)"
        )

    parser.add_argument(
        "--min_reads_per_barcode", 
        type=float, 
        help="Minimum reads per barcode (default: 1.25)", 
        default=1.25
        )
    parser.add_argument(
        "--min_covered_amplicons_per_barcode", 
        type=float, 
        help="Minimum covered amplicons per barcode (default: 0.40)", 
        default=0.40
        )
    parser.add_argument(
        "--filter_amplicons", 
        type=bool, 
        help="Filter amplicons (default: False)", 
        default=False
        )
    parser.add_argument(
        "--amplicon_variance_threshold", 
        type=float, 
        help="Amplicon variance threshold (default: 0.95)", 
        default=0.95
        )
    parser.add_argument(
        "--log_scaling", 
        type=bool, 
        help="Log scaling (default: True)", 
        default=True
        )
    parser.add_argument(
        "--thread", 
        type=int, 
        help="Number of threads (default: 1)", 
        default=1
        )
    parser.add_argument(
        "--thresholding_method", 
        type=str, 
        help="Thresholding method (default: best, options: best, kde, ellbow, completeness)", 
        default="best",
        choices=["best", "kde", "ellbow", "completeness"]
        )
    parser.add_argument(
        "--completeness_threshold", 
        type=float, 
        help="Completeness threshold (default: 0.8)", 
        default=0.8
        )
    parser.add_argument(
        "--min_cell_barcodes_count", 
        type=int, 
        help="Minimum cell barcodes count (default: 1000)", 
        default=1000
        )
    return parser.parse_args()


def main():
    args = parse_args()
    print("==>Initializing CellSavior...")
    cell_savior = CellSavior(
        args.input_barcode_distribution,
        args.input_bam,
        args.output_barcode_distribution,
        args.output_bam,
        args.output_prefix_for_plots,
        args.min_reads_per_barcode,
        args.min_covered_amplicons_per_barcode,
        args.filter_amplicons,
        args.amplicon_variance_threshold,
        args.log_scaling,
        args.thread,
        args.thresholding_method,
        args.completeness_threshold,
        args.min_cell_barcodes_count
    )
    print("==>CellSavior initialized successfully.")
    print("\n\n\n\n\n")

    print("==>Executing CellSavior...")
    cell_savior.execute()
    print("==>CellSavior executed successfully.")
    print("\n\n\n\n\n")
    
    print("==>CellSavior completed successfully.")

if __name__ == "__main__":
    main()
#--------------------------------#