<script lang="ts" module>
  // NOTE: need global cache so that drag-n-drop doesn't cause refetching (looks dumb probably)
  // TODO: make type lul
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

  import { self } from "svelte/legacy";

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

  // function pick(): void {
  //   if (bigViewCharts[0] === null) bigViewCharts[0] = myDataset;
  //   else {
  //     if (bigViewCharts[0] === myDataset) bigViewCharts[0] = null;
  //     else bigViewCharts[1] = myDataset;
  //   }
  // }
</script>

<div class="column">
  <h3 class="column-title">
    {makeTitle(dataset)}
  </h3>
  <!-- <button onclick={pick}>
    <h3
      class="dataset"
      class:selecting={bigView}
      class:selected={bigView && bigViewCharts.includes(myDataset)}
    >
      {makeTitle(dataset)}
    </h3>
  </button> -->
  {#await plots}
    <CircularProgress indeterminate style="width: 100%; height: 100%;" />
  {:then plots}
    <button onclick={() => (bigImg = plots.clustering)}>
      <img class="smol" src={plots.clustering} alt="{dataset} clustering" />
    </button>
    <button onclick={() => (bigImg = plots.violin)}>
      <img class="smol" src={plots.violin} alt="{dataset} violin" />
    </button>
  {:catch err}
    <h1>Error: {err}</h1>
  {/await}

  <!-- putting outside div breaks drag-n-drop
  <dialog open={bigImg !== null} onclick={self(() => (bigImg = null))}>
    <article class="big-modal">
      <img src={bigImg} alt="Enlarged" />
    </article>
  </dialog> -->
</div>

<style lang="scss">
  .column {
    width: 40vh;
    padding: 0.5em;
    margin: 0.5em;
    float: left;
    border: 1px solid #333333;
    overflow: hidden;
  }

  .column-title {
    margin-bottom: 1em;
    display: flex;
    justify-content: center;
    align-items: center;
  }

  img.smol {
    width: 100%;
    cursor: pointer;
  }

  // .dataset {
  //   width: 100%;
  //   text-align: center;
  //   position: relative;
  //   overflow: hidden;
  //   padding: 2px;
  //   background: var(--background-color);

  //   &.selecting {
  //     cursor: pointer;
  //     color: goldenrod;
  //     z-index: 0;

  //     &:hover,
  //     &.selected {
  //       &::before {
  //         content: "";
  //         display: block;
  //         width: calc(100% * 1.414); // sqrt(2) because diagonal
  //         padding-bottom: calc(100% * 1.414); // sqrt(2) because diagonal
  //         position: absolute;
  //         left: 50%;
  //         top: 50%;
  //         transform: translate(-50%, -50%);
  //         z-index: -2;

  //         animation: spin 8s linear infinite;
  //         background: repeating-conic-gradient(
  //           from 0deg,
  //           goldenrod 0deg 90deg,
  //           var(--background-color) 90deg 360deg
  //         );
  //       }

  //       &::after {
  //         content: "";
  //         position: absolute;
  //         inset: 2px;
  //         background: var(--background-color);
  //         z-index: -1;
  //       }
  //     }
  //   }
  // }

  // .big-modal {
  //   padding: 0;
  //   max-width: fit-content;

  //   img {
  //     height: 90vh;
  //   }
  // }
</style>
