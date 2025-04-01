<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import GeneControlsDrawerToggle from "$lib/components/GeneControls/GeneControlsDrawerToggle.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import Button, { Label } from "@smui/button";
  import { AppContent } from "@smui/drawer";
  import debounce from "lodash.debounce";
  import { onMount } from "svelte";
  import Tab from "@smui/tab";
  import TabBar from "@smui/tab-bar";

  const { token }: IndexProps = $props();

  let datasets: Dataset[] = $state([]),
    selectedDatasets: string[] = $state([]),
    loadedDatasets: boolean = $state(true),
    genes: string[] = $state([]),
    selectedGene: string = $state(""),
    loadedGenes: boolean = $state(true),
    groupBy: string = $state("Genotype"),
    splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype"),
    geneControlsOpen = $state(true),
    cellType = $state("All"),
    cellTypes = $derived(new Set(datasets.map((ds) => ds.cellType))),
    cellDatasets = $derived(
      cellType === "All"
        ? datasets
        : datasets.filter((ds) => ds.cellType === cellType)
    );
  $effect(() => {
    if (selectedDatasets.length) debounce(fetchGenes, 1000)();
    else genes = [];
  });

  function fetchGenes() {
    if (selectedDatasets.length) {
      loadedGenes = false;
      fetch(`/genes?datasets=${selectedDatasets.join(",")}`)
        .then((res) => res.json())
        .then((data) => {
          genes = data;
          loadedGenes = true;
        });
      fetch(`/preload?token=${token}&datasets=${selectedDatasets.join(",")}`, {
        method: "POST",
      });
    }
  }

  onMount(() => {
    loadedDatasets = false;
    fetch("/datasets")
      .then((res) => res.json())
      .then((data) => {
        datasets = data;
        loadedDatasets = true;
      });
  });

  const datasetColumns = [
      { key: "name", label: "Dataset" },
      { key: "year", label: "Year" },
      { key: "region", label: "Region" },
      { key: "PMID", label: "PMID" },
      { key: "species", label: "Species" },
      { key: "author", label: "Author" },
      { key: "disease", label: "Disease" },
    ],
    geneColumns = [{ key: "", label: "Gene" }];
</script>

<svelte:head>
  <title>MultiSC-Viewer</title>
</svelte:head>

<Navbar />
<div style="display: inline-flex; gap: 1rem;">
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
