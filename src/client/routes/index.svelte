<script lang="ts">
  import DataTableSearch from "$lib/components/DataTableSearch.svelte";
  import GeneControls from "$lib/components/GeneControls.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import Button, { Label } from "@smui/button";
  import debounce from "lodash.debounce";
  import { onMount } from "svelte";

  const { token }: IndexProps = $props();

  let datasets: object[] = $state([]),
    selectedDatasets: string[] = $state([]),
    loadedDatasets: boolean = $state(true),
    genes: object[] = $state([]),
    selectedGene: string = $state(""),
    loadedGenes: boolean = $state(true),
    groupBy: string = $state("Genotype"),
    splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype");

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
<div style="display: flex;align-items: center;">
  <GeneControls bind:groupBy />
  <Button onclick={plot} variant="raised" disabled={!selectedDatasets.length}>
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
<DataTableSearch
  label="Genes"
  data={genes}
  columns={geneColumns}
  bind:selected={selectedGene}
  loaded={loadedGenes}
/>
