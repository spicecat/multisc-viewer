<script lang="ts">
  import DataTable from "$lib/components/DataTable.svelte";
  import IconButton from "@smui/icon-button";
  import TextField from "@smui/textfield";
  import { Cell, InnerGrid } from "@smui/layout-grid";
  import Fuse from "fuse.js";

  let {
    label = "",
    data = [],
    columns = [],
    selected = $bindable(),
    searchOptions = {},
    loaded = true,
  } = $props();

  let defaultSearchOptions = $derived({
      keys: columns.map(({ key }) => key),
      threshold: 0.0,
      ignoreLocation: true,
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

<div style="display: inline-block;">
  <div>
    <TextField bind:value={query} {label}>
      {#snippet leadingIcon()}
        <IconButton class="material-icons">search</IconButton>
      {/snippet}
    </TextField>
  </div>
  <DataTable data={filteredItems} {columns} bind:selected {loaded} />
</div>
