<script lang="ts">
  import DatasetsTable from "$lib/components/DataTable/DatasetsTable.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Button, { Label } from "@smui/button";
  import { AppContent } from "@smui/drawer";
  import debounce from "lodash.debounce";

  let { datasets }: { datasets: Dataset[] } = $props();

  let selectedDatasets: string[] = $state([]);
  let genes: string[] = $state([]);
  let isLoadingGenes: boolean = $state(true);
  let selectedGene: string = $state("");
  let groupBy: Grouping = $state("Genotype");
  let geneControlsOpen: boolean = $state(true);

  const splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype");

  $effect(() => {
    if (selectedDatasets.length) debounce(fetchGenes, 1000)();
    else genes = [];
  });

  async function fetchGenes() {
    if (selectedDatasets.length == 0) return;

    isLoadingGenes = true;
    try {
      const response = await fetch(
        `/genes?datasets=${selectedDatasets.join(",")}`
      );
      genes = await response.json();
    } catch (error) {
      console.error("Error fetching genes:", error);
    } finally {
      isLoadingGenes = false;
    }
  }
</script>

<div style="display: flex;">
  <GeneControlsDrawer
    {genes}
    isLoading={isLoadingGenes}
    bind:selected={selectedGene}
    bind:groupBy
    {geneControlsOpen}
  />

  <AppContent style="margin-left: 0;">
    <div style="align-items: center;display: flex;gap: 1rem;">
      <GeneControlsDrawerToggle bind:geneControlsOpen />
      {#if selectedDatasets.length}
        <Button
          variant="raised"
          href={`/compare?datasets=${selectedDatasets.join(",")}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`}
        >
          <Label>Plot</Label>
        </Button>
      {:else}
        <Button variant="raised" {...{ disabled: true }}>
          <Label>Plot</Label>
        </Button>
      {/if}
    </div>

    <DatasetsTable {datasets} bind:selected={selectedDatasets} />
  </AppContent>
</div>
