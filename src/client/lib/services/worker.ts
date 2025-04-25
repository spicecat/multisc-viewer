import { EventEmitter } from '../utils/events';
import type { WorkerMessage, DataFilter } from '../types';

class WorkerService extends EventEmitter {
  private worker: Worker;
  private pendingTasks = new Map<string, {
    resolve: (value: any) => void;
    reject: (error: any) => void;
  }>();

  constructor() {
    super();
    this.worker = new Worker(new URL('../workers/main.ts', import.meta.url));
    this.setupWorker();
  }

  private setupWorker() {
    this.worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
      const { id, type, payload } = event.data;

      if (type === 'error') {
        this.handleError(id, payload);
      } else {
        this.handleSuccess(id, payload);
      }
    };

    this.worker.onerror = (error) => {
      console.error('Worker error:', error);
      this.emit('worker:error', error);
    };
  }

  private handleSuccess(id: string, payload: any) {
    const task = this.pendingTasks.get(id);
    if (task) {
      task.resolve(payload);
      this.pendingTasks.delete(id);
    }
  }

  private handleError(id: string, error: any) {
    const task = this.pendingTasks.get(id);
    if (task) {
      task.reject(new Error(error.message || 'Worker task failed'));
      this.pendingTasks.delete(id);
    }
  }

  private generateTaskId(): string {
    return Math.random().toString(36).substring(2, 15);
  }

  async executeTask<T>(type: string, payload: any): Promise<T> {
    const id = this.generateTaskId();

    return new Promise((resolve, reject) => {
      this.pendingTasks.set(id, { resolve, reject });

      this.worker.postMessage({
        id,
        type,
        payload
      } as WorkerMessage);
    });
  }

  // Pre-defined task types
  async processData<T>(data: T[]): Promise<T[]> {
    return this.executeTask('processData', data);
  }

  async calculateStatistics<T>(data: T[]): Promise<{
    count: number;
    numeric: Record<string, {
      min: number;
      max: number;
      mean: number;
      median: number;
    }>;
    categorical: Record<string, {
      categories: string[];
      counts: Record<string, number>;
    }>;
  }> {
    return this.executeTask('calculateStatistics', data);
  }

  async filterDataset<T>(data: T[], filters: DataFilter): Promise<T[]> {
    return this.executeTask('filterDataset', { data, filters });
  }

  terminate(): void {
    this.worker.terminate();
    this.pendingTasks.clear();
    this.emit('terminated');
  }
}

// Export singleton instance
export const workerService = new WorkerService();