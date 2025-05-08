<script lang="ts">
  import DataTableSearch from "./DataTableSearch.svelte";

  let {
    genes,
    isLoading,
    selected = $bindable(),
  }: { genes: Gene[]; isLoading: boolean; selected: string } = $props();

  const geneColumns: Column[] = [{ key: "", label: "Gene" }];
</script>

<DataTableSearch
  label="Gene"
  data={genes}
  columns={geneColumns}
  {isLoading}
  bind:selected
/>
