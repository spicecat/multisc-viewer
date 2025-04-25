import { api } from './api';
import { workerService } from './worker';
import { createAsyncStore } from '../utils/store';
import { showLoading, hideLoading } from '../stores/uiState';
import type { Publication, Dataset } from '../../../interfaces/types';
import type { PublicationMetrics, DataFilter } from '../types';

class PublicationService {
  private readonly publicationStore = createAsyncStore<Publication[]>();
  private readonly metricsCache = new Map<string, PublicationMetrics>();

  async loadPublications(): Promise<Publication[]> {
    showLoading('Loading publications...');
    try {
      const publications = await api.getPublications();
      this.publicationStore.load(Promise.resolve(publications));
      hideLoading();
      return publications;
    } catch (error) {
      hideLoading();
      throw error;
    }
  }

  async getPublicationMetrics(publication: Publication, datasets: Dataset[]): Promise<PublicationMetrics> {
    const cacheKey = publication.publicationId;
    if (this.metricsCache.has(cacheKey)) {
      return this.metricsCache.get(cacheKey)!;
    }

    const publicationDatasets = datasets.filter(d => d.publication?.id === publication.publicationId);
    
    const metrics: PublicationMetrics = {
      datasetCount: publicationDatasets.length,
      totalCells: publicationDatasets.reduce((sum, d) => sum + (d.cells?.length ?? 0), 0),
      totalGenes: publicationDatasets.reduce((sum, d) => sum + (d.genes?.length ?? 0), 0),
      uniqueGenotypes: new Set(
        publicationDatasets.flatMap(d => d.genotypes || [])
      ),
      uniqueCellTypes: new Set(
        publicationDatasets.flatMap(d => d.cellTypes || [])
      )
    };

    this.metricsCache.set(cacheKey, metrics);
    return metrics;
  }

  async analyzePublication(publication: Publication, datasets: Dataset[]): Promise<any> {
    const publicationDatasets = datasets.filter(d => d.publication?.id === publication.publicationId);
    const combinedData = publicationDatasets.flatMap(d => d.data || []);
    
    return workerService.calculateStatistics(combinedData);
  }

  async comparePublications(publications: Publication[], datasets: Dataset[]): Promise<any> {
    const publicationComparisons = await Promise.all(
      publications.map(async publication => {
        const metrics = await this.getPublicationMetrics(publication, datasets);
        const analysis = await this.analyzePublication(publication, datasets);
        
        return {
          publication,
          metrics,
          analysis
        };
      })
    );

    return workerService.processData(publicationComparisons);
  }

  async filterPublications(
    publications: Publication[],
    filters: DataFilter
  ): Promise<Publication[]> {
    return workerService.filterDataset(publications, filters);
  }

  subscribe(callback: (state: { data: Publication[] | null }) => void) {
    return this.publicationStore.data.subscribe(callback);
  }

  clearCache(): void {
    this.metricsCache.clear();
  }
}

// Export singleton instance
export const publicationService = new PublicationService();