<script lang="ts">
  import DataTable from "$lib/components/DataTable.svelte";
  import IconButton from "@smui/icon-button";
  import TextField from "@smui/textfield";
  import Fuse from "fuse.js";

  let {
    label = "",
    data = [],
    columns = [],
    selected = $bindable(),
    searchOptions = {},
  } = $props();

  let defaultSearchOptions = $derived({
      keys: columns.map(({ key }) => key),
      useExtendedSearch: true,
    }),
    fuse = $derived(
      new Fuse(data, {
        ...defaultSearchOptions,
        ...searchOptions,
      })
    ),
    query: string = $state(""),
    filteredItems = $derived(
      query ? fuse.search(query).map(({ item }) => item) : data
    );
</script>

<TextField bind:value={query} {label}>
  {#snippet leadingIcon()}
    <IconButton class="material-icons">search</IconButton>
  {/snippet}
</TextField>
<DataTable data={filteredItems} {columns} bind:selected />
