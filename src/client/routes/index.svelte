<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import GeneGroupSplit from "$lib/components/GeneControls/GeneGroupSplit.svelte";
  import GeneControlsDrawer from "$lib/components/GeneControls/GeneControlsDrawer.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import Button, { Label } from "@smui/button";
  import debounce from "lodash.debounce";
  import { onMount } from "svelte";
  import Drawer, {
    AppContent,
    Content,
    Header,
    Title,
    Scrim,
  } from "@smui/drawer";

  const { token }: IndexProps = $props();

  let datasets: object[] = $state([]),
    selectedDatasets: string[] = $state([]),
    loadedDatasets: boolean = $state(true),
    genes: object[] = $state([]),
    selectedGene: string = $state(""),
    loadedGenes: boolean = $state(true),
    groupBy: string = $state("Genotype"),
    splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype"),
    geneControlsOpen = $state(true);
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

  function plot(): void {
    if (selectedDatasets.length)
      window.location.href = `/compare?datasets=${selectedDatasets.join(",")}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`; // use goto
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
      { key: "name", label: "Name" },
      { key: "year", label: "Year", sortable: true },
      { key: "region", label: "Region" },
      { key: "PMID", label: "PMID" },
      { key: "species", label: "Species" },
      { key: "author", label: "Author" },
      { key: "disease", label: "Disease" },
    ],
    geneColumns = [{ key: "", label: "Gene" }];
</script>

<svelte:head>
  <title>Dataset Comparison</title>
</svelte:head>

<Navbar />
<div style="display: flex;">
  <GeneControlsDrawer
    {genes}
    bind:selectedGene
    {loadedGenes}
    bind:groupBy
    {geneControlsOpen}
  />
  <AppContent>
    <!-- <GeneGroupSplit bind:groupBy /> -->
    <div>
      <Button
        onclick={() => (geneControlsOpen = !geneControlsOpen)}
        variant="raised"
      >
        <Label>Open Gene Controls</Label>
      </Button>
      <Button
        variant="raised"
        {...selectedDatasets.length
          ? {
              href: `/compare?datasets=${selectedDatasets.join(",")}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`,
            }
          : { disabled: true }}
      >
        <Label>Plot</Label>
      </Button>
      <Button
        variant="raised"
        href="https://google.com"
        {...selectedDatasets.length
          ? {
              href: `/compare?datasets=${selectedDatasets.join(",")}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`,
            }
          : {
            href: "javascript:void(0);",
            disabled: true
             }}
      >
        <Label>Plot</Label>
      </Button>
    </div>
    <DataTableSearch
      label="Datasets"
      data={datasets}
      columns={datasetColumns}
      bind:selected={selectedDatasets}
      loaded={loadedDatasets}
    />
    <!-- <DataTableSearch
    label="Genes"
    data={genes}
    columns={geneColumns}
    bind:selected={selectedGene}
    loaded={loadedGenes}
    /> -->
  </AppContent>
</div>
