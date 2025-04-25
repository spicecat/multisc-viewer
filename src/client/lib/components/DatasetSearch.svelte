<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Button, { Label } from "@smui/button";
  import { AppContent } from "@smui/drawer";
  import Tab from "@smui/tab";
  import TabBar from "@smui/tab-bar";
  import debounce from "lodash.debounce";
  import type { Dataset } from "../../interfaces/types";

  // Props
  const { datasets = [] } = $props();

  // State
  let selectedDatasets: string[] = $state([]);
  let loadedDatasets: boolean = $state(true);
  let genes: string[] = $state([]);
  let selectedGene: string = $state("");
  let loadedGenes: boolean = $state(true);
  let groupBy: string = $state("Genotype");
  let geneControlsOpen = $state(true);
  let cellType = $state("All");

  // Derived state
  const splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype");
  const cellTypes = $derived(new Set(datasets.map((ds) => ds.cellType)));
  const cellDatasets = $derived(
    cellType === "All"
      ? datasets
      : datasets.filter((ds) => ds.cellType === cellType)
  );

  // Dataset table configuration
  const datasetColumns = [
    { key: "name", label: "Dataset" },
    { key: "year", label: "Year" },
    { key: "region", label: "Region" },
    { key: "PMID", label: "PMID" },
    { key: "species", label: "Species" },
    { key: "author", label: "Author" },
    { key: "disease", label: "Disease" },
    { key: "cellType", label: "Cell Type" },
  ];

  // Effects and handlers
  $effect(() => {
    if (selectedDatasets.length) {
      debounce(fetchGenes, 1000)();
    } else {
      genes = [];
    }
  });

  async function fetchGenes() {
    if (!selectedDatasets.length) return;

    loadedGenes = false;
    try {
      const response = await fetch(
        `/genes?datasets=${selectedDatasets.join(",")}`
      );
      genes = await response.json();
    } catch (error) {
      console.error("Error fetching genes:", error);
    } finally {
      loadedGenes = true;
    }
  }
</script>

<div>
  <GeneControlsDrawer
    {genes}
    bind:selectedGene
    {loadedGenes}
    bind:groupBy
    {geneControlsOpen}
  />

  <AppContent>
    <div>
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

    <TabBar tabs={["All", ...cellTypes]} bind:active={cellType}>
      {#snippet tab(tab)}
        <Tab {tab}>
          <Label>{tab}</Label>
        </Tab>
      {/snippet}
    </TabBar>

    <DataTableSearch
      label="Datasets"
      data={cellDatasets}
      columns={datasetColumns}
      bind:selected={selectedDatasets}
      loaded={loadedDatasets}
    />
  </AppContent>
</div>
