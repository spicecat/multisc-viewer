<script lang="ts">
  import ChartDisplay from "$lib/components/ChartDisplay.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import IconButton from "@smui/icon-button";
  import Drawer, {
    AppContent,
    Content,
    Header,
    Title,
    Scrim,
  } from "@smui/drawer";
  import meta from "$meta";
  import Button, { Label } from "@smui/button";
  import { dndzone } from "svelte-dnd-action";
  import { self } from "svelte/legacy";

  const { genes, gene, order }: CompareProps = $props();
  const { query } = meta;
  console.log(genes, gene, order);

  let groupBy: string = $state("Genotype"),
    selectedGene: string = $state(gene),
    dsOrder = $state(order.map((ds) => ({ id: ds }))),
    bigView: boolean = $state(false),
    bigViewCharts: (Promise<RenderResult> | null)[] = $state([null, null]),
    geneControlsOpen = $state(false);

  const splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype"),
    config = $derived({ selectedGene, groupBy, splitBy });

  $effect(() => {
    if (
      !history.state ||
      history.state.selectedGene !== selectedGene ||
      history.state.groupBy !== groupBy
    )
      history.pushState(
        { selectedGene, groupBy },
        "",
        `/compare?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`
      );
  });
  const geneColumns = [{ key: "", label: "Gene" }];
</script>

<Navbar />

<div style="display: flex;">
  <GeneControlsDrawer
    {genes}
    bind:selectedGene
    loadedGenes={true}
    bind:groupBy
    {geneControlsOpen}
  />
  <Scrim />
  <AppContent>
    <div>
      <GeneControlsDrawerToggle bind:geneControlsOpen />
    </div>
    <div
      class="row charts"
      use:dndzone={{ items: dsOrder }}
      onconsider={(evt) => (dsOrder = evt.detail.items)}
      onfinalize={(evt) => (dsOrder = evt.detail.items)}
    >
      {#each dsOrder as { id: dataset } (dataset)}
        <ChartDisplay {dataset} {config} {bigView} {bigViewCharts} />
      {/each}
    </div>
  </AppContent>
</div>

<!-- <dialog
  open={!!(bigView && bigViewCharts[0] && bigViewCharts[1])}
  onclick={self(() => {
    bigView = false;
    bigViewCharts = [null, null];
  })}
>
  <article class="big-modal">
    {#if bigView && bigViewCharts[0] && bigViewCharts[1]}
      <div class="row">
        <div class="col center">
          {#await bigViewCharts[0]}
            Loading...
          {:then charts}
            <img src={charts.clustering} alt="Enlarged" />
            <img src={charts.violin} alt="Enlarged" />
          {/await}
        </div>
        <div class="col center">
          {#await bigViewCharts[1]}
            Loading...
          {:then charts}
            <img src={charts.clustering} alt="Enlarged" />
            <img src={charts.violin} alt="Enlarged" />
          {/await}
        </div>
      </div>
    {/if}
  </article>
</dialog> -->

<style lang="scss">
  @keyframes spin {
    0% {
      transform: rotate(0deg);
    }

    100% {
      transform: rotate(360deg);
    }
  }

  :global(.fa-spinner) {
    animation: spin 1s linear infinite;
  }

  .charts {
    max-width: 100vw;
    display: flex;
  }

  .big-modal {
    padding: 0;
    max-width: fit-content;

    .row {
      gap: 5vh;

      .col {
        gap: 5vh;
        width: 40vw;

        img {
          height: 40vh;
        }
      }
    }
  }

</style>
