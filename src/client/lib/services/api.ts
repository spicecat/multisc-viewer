import { handleHttpError } from '../utils/errors';
import type { Publication, Dataset } from '../../../interfaces/types';
import type { ApiOptions } from '../types';

class ApiService {
  private baseUrl: string;
  private headers: Record<string, string>;

  constructor(options: ApiOptions = {}) {
    this.baseUrl = options.baseUrl || '';
    this.headers = {
      'Content-Type': 'application/json',
      ...options.headers
    };
  }

  private async request<T>(
    endpoint: string,
    options: RequestInit = {}
  ): Promise<T> {
    const response = await fetch(`${this.baseUrl}${endpoint}`, {
      ...options,
      headers: {
        ...this.headers,
        ...options.headers
      }
    });

    if (!response.ok) {
      handleHttpError(response);
    }

    return response.json();
  }

  // Publication endpoints
  async getPublications(): Promise<Publication[]> {
    return this.request('/publication');
  }

  // Dataset endpoints
  async getDatasets(): Promise<Dataset[]> {
    return this.request('/datasets');
  }

  async getGenes(datasets: string[]): Promise<string[]> {
    return this.request(`/genes?datasets=${datasets.join(',')}`);
  }

  async renderPlot(
    dataset: string,
    gene: string,
    groupBy: string,
    splitBy: string,
    token: string | null = null
  ): Promise<{ clustering: string; violin: string }> {
    const params = new URLSearchParams({
      dataset,
      gene,
      groupBy,
      splitBy,
      ...(token && { token })
    });

    return this.request(`/plot?${params}`);
  }

  async preloadDatasets(
    datasets: string[],
    token: string | null = null
  ): Promise<void> {
    const params = new URLSearchParams({
      datasets: datasets.join(','),
      ...(token && { token })
    });

    return this.request('/preload', {
      method: 'POST',
      body: JSON.stringify({ datasets }),
    });
  }
}

// Create and export a singleton instance
export const api = new ApiService();