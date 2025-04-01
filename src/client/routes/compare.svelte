<script lang="ts">
  import ChartDisplay from "$lib/components/ChartDisplay.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import meta from "$meta";
  import { AppContent, Scrim } from "@smui/drawer";
  import { dndzone } from "svelte-dnd-action";

  const { genes, gene, order }: CompareProps = $props();
  const { query } = meta;

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
</script>

<svelte:head>
  <title>MultiSC-Viewer Plot</title>
</svelte:head>

<Navbar />
<div style="display: inline-flex; gap: 1rem;">
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
    <section
      class="board"
      use:dndzone={{ items: dsOrder }}
      onconsider={(evt) => (dsOrder = evt.detail.items)}
      onfinalize={(evt) => (dsOrder = evt.detail.items)}
    >
      {#each dsOrder as { id: dataset } (dataset)}
        <ChartDisplay {dataset} {config} {bigView} {bigViewCharts} />
      {/each}
    </section>
  </AppContent>
</div>

<style lang="scss">
  .board {
    min-height: 40vh;
    margin: 0.5em;
    display: flex;
    overflow-x: scroll;
  }

  // .big-modal {
  //   padding: 0;
  //   max-width: fit-content;

  //   .row {
  //     gap: 5vh;

  //     .col {
  //       gap: 5vh;
  //       width: 40vw;

  //       img {
  //         height: 40vh;
  //       }
  //     }
  //   }
  // }
</style>
