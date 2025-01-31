<script lang="ts">
  import { self } from "svelte/legacy";

  const { dataset, config }: { dataset: string; config: PlotConfig } = $props();

  const { selectedGene, groupBy, splitBy } = $derived(config);

  let bigImg: string | null = $state(null);
</script>

<div class="col">
  <h3 class="dataset">{dataset.replaceAll("_", " ")}</h3>
  {#await fetch(`/plot?dataset=${dataset}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`).then( (res) => res.json() )}
    <div>
      <i class="fa-solid fa-spinner"></i> Generating Plots...
    </div>
  {:then plots}
    <img
      class="smol"
      src={plots.clustering}
      alt="{dataset} clustering"
      onclick={() => (bigImg = plots.clustering)}
    />
    <img
      class="smol"
      src={plots.violin}
      alt="{dataset} violin"
      onclick={() => (bigImg = plots.violin)}
    />
  {:catch err}
    <h1>Error: {err}</h1>
  {/await}
</div>

<dialog open={bigImg !== null} onclick={self(() => (bigImg = null))}>
  <article class="big-modal">
    <img src={bigImg} alt="Enlarged" />
  </article>
</dialog>

<style lang="scss">
  .col,
  img.smol {
    width: 300px;
  }

  .col {
    flex-shrink: 0;
  }

  img.smol {
    cursor: pointer;
  }

  .dataset {
    text-align: center;
  }

  .big-modal {
    padding: 0;
    max-width: fit-content;

    img {
      height: 90vh;
    }
  }
</style>
