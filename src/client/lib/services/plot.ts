import { api } from './api';
import { createStore } from '../utils/store';
import { showLoading, hideLoading } from '../stores/uiState';
import { appConfig } from '../../../config/app.config';
import type { PlotOptions } from '../types';

interface PlotCache {
  [key: string]: {
    data: { clustering: string; violin: string };
    timestamp: number;
  };
}

interface PlotRequest {
  dataset: string;
  gene: string;
  groupBy: string;
  splitBy: string;
  options?: PlotOptions;
}

class PlotService {
  private static MAX_CACHE_SIZE = 50; // Maximum number of plots to cache
  private static CACHE_TTL = 1000 * 60 * 30; // 30 minutes
  private static DEBOUNCE_DELAY = 250; // ms

  private cache: PlotCache = {};
  private cacheQueue: string[] = []; // LRU tracking
  private pendingRequests = new Map<string, Promise<any>>();
  private debounceTimers = new Map<string, NodeJS.Timeout>();

  private readonly plotStore = createStore<{
    plots: Map<string, { clustering: string; violin: string }>;
  }>({
    plots: new Map()
  });

  private getCacheKey(request: PlotRequest): string {
    const { dataset, gene, groupBy, splitBy, options = {} } = request;
    const optionsKey = Object.entries(options)
      .sort(([a], [b]) => a.localeCompare(b))
      .map(([k, v]) => `${k}:${JSON.stringify(v)}`)
      .join(',');
    
    return `${dataset}-${gene}-${groupBy}-${splitBy}${optionsKey ? `-${optionsKey}` : ''}`;
  }

  private isCacheValid(timestamp: number): boolean {
    return Date.now() - timestamp < PlotService.CACHE_TTL;
  }

  private evictStaleCache() {
    const now = Date.now();
    Object.entries(this.cache).forEach(([key, value]) => {
      if (!this.isCacheValid(value.timestamp)) {
        delete this.cache[key];
        this.cacheQueue = this.cacheQueue.filter(k => k !== key);
      }
    });
  }

  private addToCache(key: string, data: any) {
    // Evict old entries if we're at capacity
    while (this.cacheQueue.length >= PlotService.MAX_CACHE_SIZE) {
      const oldestKey = this.cacheQueue.shift();
      if (oldestKey) delete this.cache[oldestKey];
    }

    this.cache[key] = {
      data,
      timestamp: Date.now()
    };
    this.cacheQueue.push(key);
  }

  async renderPlot(request: PlotRequest): Promise<{ clustering: string; violin: string }> {
    const cacheKey = this.getCacheKey(request);

    // Cancel any pending debounced requests for this key
    if (this.debounceTimers.has(cacheKey)) {
      clearTimeout(this.debounceTimers.get(cacheKey));
    }

    // Debounce new requests
    return new Promise((resolve, reject) => {
      const timer = setTimeout(async () => {
        try {
          const result = await this.executeRenderPlot(request);
          resolve(result);
        } catch (error) {
          reject(error);
        }
      }, PlotService.DEBOUNCE_DELAY);
      this.debounceTimers.set(cacheKey, timer);
    });
  }

  private async executeRenderPlot(request: PlotRequest, retryCount = 0): Promise<{ clustering: string; violin: string }> {
    const cacheKey = this.getCacheKey(request);
    
    // Evict stale cache entries
    this.evictStaleCache();

    // Check cache first
    const cached = this.cache[cacheKey];
    if (cached && this.isCacheValid(cached.timestamp)) {
      // Move to end of LRU queue
      this.cacheQueue = this.cacheQueue.filter(k => k !== cacheKey);
      this.cacheQueue.push(cacheKey);
      return cached.data;
    }

    // Check if there's a pending request
    if (this.pendingRequests.has(cacheKey)) {
      return this.pendingRequests.get(cacheKey)!;
    }

    // Show loading state
    showLoading('Rendering plot...');

    // Create new request with retry logic
    const requestPromise = api.renderPlot(
      request.dataset,
      request.gene,
      request.groupBy,
      request.splitBy
    )
    .then(result => {
      // Update cache with LRU
      this.addToCache(cacheKey, result);

      // Update store
      this.plotStore.update(state => {
        const plots = new Map(state.plots);
        plots.set(cacheKey, result);
        return { plots };
      });

      // Cleanup
      this.pendingRequests.delete(cacheKey);
      hideLoading();

      return result;
    })
    .catch(async error => {
      this.pendingRequests.delete(cacheKey);
      hideLoading();

      // Retry logic for recoverable errors
      if (retryCount < 3 && error.status !== 404) {
        await new Promise(resolve => setTimeout(resolve, Math.pow(2, retryCount) * 1000));
        return this.executeRenderPlot(request, retryCount + 1);
      }
      
      throw error;
    });

    // Store pending request
    this.pendingRequests.set(cacheKey, requestPromise);

    return requestPromise;
  }

  clearCache(): void {
    this.cache = {};
    this.cacheQueue = [];
    this.plotStore.update(() => ({ plots: new Map() }));
  }

  subscribe(callback: (state: { plots: Map<string, { clustering: string; violin: string }> }) => void) {
    return this.plotStore.subscribe(callback);
  }
}

// Export singleton instance
export const plotService = new PlotService();