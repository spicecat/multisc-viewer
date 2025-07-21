<script lang="ts">
  import Tab, { Label } from "@smui/tab";
  import TabBar from "@smui/tab-bar";
  import DataTableSearch from "./DataTableSearch.svelte";

  let {
    datasets,
    selected = $bindable(),
  }: {
    datasets: Dataset[];
    selected: string[];
  } = $props();

  let cellType = $state("");

  const cellTypes = $derived([...new Set(datasets.map((ds) => ds.cellType))]);
  const cellDatasets = $derived(
    datasets.filter((ds) => ds.cellType === cellType)
  );

  $effect(() => {
    if (cellTypes.length > 0 && !cellType) cellType = cellTypes[0];
  });

  const datasetColumns: Column[] = [
    { key: "name", label: "Dataset" },
    { key: "year", label: "Year" },
    { key: "region", label: "Region" },
    { key: "PMID", label: "PMID" },
    { key: "species", label: "Species" },
    { key: "author", label: "Author" },
    { key: "disease", label: "Disease" },
    { key: "cellType", label: "Cell Type" },
    // { key: "publication", label: "Publication", url: true }, // TODO
  ];
</script>

<TabBar tabs={cellTypes} bind:active={cellType}>
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
  bind:selected
/>
