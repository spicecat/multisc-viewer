<script lang="ts" module>
  const cache = new Map<string, Promise<RenderResult>>();

  function cachedFetch(url: string): Promise<RenderResult> {
    if (cache.has(url)) return cache.get(url);
    else {
      const res = fetch(url).then((res) => res.json());
      cache.set(url, res);
      return res;
    }
  }
</script>

<script lang="ts">
  import { makeTitle } from "$lib/utils/utils";
  import CircularProgress from "@smui/circular-progress";

  const {
    dataset,
    config,
    bigView,
    bigViewCharts,
  }: {
    dataset: string;
    config: PlotConfig;
    bigView: boolean;
    bigViewCharts: (Promise<RenderResult> | null)[];
  } = $props();
  const ds = dataset;

  const { selectedGene, groupBy, splitBy } = $derived(config),
    plots = $derived(
      cachedFetch(
        `/plot?dataset=${dataset}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`
      )
    );

  let bigImg: string | null = $state(null);
</script>

<div class="chart-column">
  <h3 class="chart-title">
    {makeTitle(dataset)}
  </h3>

  {#await plots}
    <div class="loading-container">
      <CircularProgress indeterminate />
    </div>
  {:then plotData}
    <div class="chart-images">
      <button
        class="chart-button"
        onclick={() => (bigImg = plotData.clustering)}
        aria-label="View clustering plot"
      >
        <img
          class="chart-image"
          src={plotData.clustering}
          alt="{dataset} clustering"
        />
      </button>

      <button
        class="chart-button"
        onclick={() => (bigImg = plotData.violin)}
        aria-label="View violin plot"
      >
        <img class="chart-image" src={plotData.violin} alt="{dataset} violin" />
      </button>
    </div>
  {:catch error}
    <div class="error-message">
      <p>Error: {error.message}</p>
    </div>
  {/await}
</div>

<style lang="scss">
  .chart-column {
    width: 40vh;
    padding: 0.8rem;
    margin: 0.5rem;
    border: 1px solid var(--border-color, #333);
    border-radius: 4px;
    background-color: var(--background-color, #fff);
  }

  .chart-title {
    margin-bottom: 1rem;
    text-align: center;
    font-size: 1.1rem;
    font-weight: 500;
  }

  .loading-container {
    display: flex;
    justify-content: center;
    align-items: center;
    min-height: 200px;
  }

  .chart-images {
    display: flex;
    flex-direction: column;
    gap: 1rem;
  }

  .chart-button {
    padding: 0;
    border: none;
    background: none;
    cursor: pointer;
    transition: transform 0.2s ease;

    &:hover {
      transform: scale(1.02);
    }

    &:focus {
      outline: 2px solid var(--primary-color, #2196f3);
      outline-offset: 2px;
    }
  }

  .chart-image {
    width: 100%;
    height: auto;
    border-radius: 4px;
  }

  .error-message {
    color: var(--error-color, #f44336);
    text-align: center;
    padding: 1rem;
  }
</style>
