<script lang="ts">
  import ChartDisplay from "$lib/components/ChartDisplay.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import meta from "$meta";
  import { AppContent, Scrim } from "@smui/drawer";
  import Button from "@smui/button";
  import { dndzone } from "svelte-dnd-action";
  import html2canvas from "html2canvas";
  import { tick } from "svelte";

  interface RenderResult {
    clustering: string;
    violin: string;
  }
  interface CompareProps {
    genes: string[];
    gene: string;
    order: string[];
  }

  const { genes, gene, order }: CompareProps = $props();
  const { query } = meta;

  let groupBy: string = $state("Genotype"),
    selectedGene: string = $state(gene),
    dsOrder = $state(order.map((ds) => ({ id: ds }))),
    bigView: boolean = $state(false),
    bigViewCharts: (Promise<RenderResult> | null)[] = $state([null, null]),
    geneControlsOpen = $state(false),
    boardElement: HTMLElement | null = $state(null),
    isDownloading = $state(false);

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

  async function downloadBoard() {
    if (!boardElement || isDownloading) return;

    isDownloading = true;
    await tick();

    try {
      const canvas = await html2canvas(boardElement, {
        useCORS: true,
        allowTaint: true,
        scrollX: 0,
        scrollY: 0,
        windowWidth: boardElement.scrollWidth,
        windowHeight: boardElement.scrollHeight,
      });

      const image = canvas.toDataURL("image/png");

      const link = document.createElement("a");
      link.href = image;
      link.download = `MultiSC-Viewer-${selectedGene}-${groupBy}.png`;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    } catch (error) {
      console.error("Error downloading board:", error);
    } finally {
      isDownloading = false;
    }
  }
</script>

<svelte:head>
  <title>MultiSC-Viewer Plot</title>
</svelte:head>

<Navbar />
<main style="display: inline-flex; gap: 1rem;">
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
      <Button variant="raised" onclick={downloadBoard}>
        {#if isDownloading}
          Downloading...
        {:else}
          Download
        {/if}
      </Button>
    </div>
    <section
      class="board"
      bind:this={boardElement}
      use:dndzone={{ items: dsOrder }}
      onconsider={(evt) => (dsOrder = evt.detail.items)}
      onfinalize={(evt) => (dsOrder = evt.detail.items)}
    >
      {#each dsOrder as { id: dataset } (dataset)}
        <ChartDisplay {dataset} {config} {bigView} {bigViewCharts} />
      {/each}
    </section>
  </AppContent>
</main>

<style lang="scss">
  .board {
    min-height: 40vh;
    margin: 0.5em;
    display: flex;
    overflow-x: scroll;
    background-color: white;
    padding: 1em;
  }
</style>
