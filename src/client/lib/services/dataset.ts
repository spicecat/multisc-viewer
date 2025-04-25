import { api } from './api';
import { workerService } from './worker';
import { createAsyncStore } from '../utils/store';
import { showLoading, hideLoading } from '../stores/uiState';
import type { Dataset } from '../../../interfaces/types';
import type { DatasetMetrics, DataFilter } from '../types';

class DatasetService {
  private readonly datasetStore = createAsyncStore<Dataset[]>();
  private readonly geneStore = createAsyncStore<Map<string, string[]>>();
  private metricsCache = new Map<string, DatasetMetrics>();

  async loadDatasets(): Promise<Dataset[]> {
    showLoading('Loading datasets...');
    try {
      const datasets = await api.getDatasets();
      this.datasetStore.load(Promise.resolve(datasets));
      hideLoading();
      return datasets;
    } catch (error) {
      hideLoading();
      throw error;
    }
  }

  async loadGenes(datasets: string[]): Promise<string[]> {
    const cacheKey = datasets.sort().join(',');
    const cached = this.geneStore.data?.get(cacheKey);
    if (cached) {
      return cached;
    }

    showLoading('Loading genes...');
    try {
      const genes = await api.getGenes(datasets);
      this.geneStore.data?.set(cacheKey, genes);
      hideLoading();
      return genes;
    } catch (error) {
      hideLoading();
      throw error;
    }
  }

  async getMetrics(dataset: Dataset): Promise<DatasetMetrics> {
    if (this.metricsCache.has(dataset.name)) {
      return this.metricsCache.get(dataset.name)!;
    }

    const metrics: DatasetMetrics = {
      geneCount: dataset.genes?.length ?? 0,
      cellCount: dataset.cells?.length ?? 0,
      genotypeCount: new Set(dataset.genotypes).size,
      cellTypeCount: new Set(dataset.cellTypes).size
    };

    this.metricsCache.set(dataset.name, metrics);
    return metrics;
  }

  async filterDatasets(
    datasets: Dataset[],
    filters: DataFilter
  ): Promise<Dataset[]> {
    return workerService.filterDataset(datasets, filters);
  }

  async analyzeDataset(dataset: Dataset): Promise<any> {
    return workerService.calculateStatistics(dataset.data || []);
  }

  subscribe(callback: (state: { data: Dataset[] | null }) => void) {
    return this.datasetStore.data.subscribe(callback);
  }

  subscribeToGenes(callback: (state: { data: Map<string, string[]> | null }) => void) {
    return this.geneStore.data.subscribe(callback);
  }
}

// Export singleton instance
export const datasetService = new DatasetService();