<script lang="ts">
  import DataTable from "./DataTable.svelte";
  import IconButton from "@smui/icon-button";
  import TextField from "@smui/textfield";
  import Fuse from "fuse.js";

  interface Column {
    key: string;
    label: string;
  }

  // Props with defaults
  let {
    label = "",
    data = [],
    columns = [] as Column[],
    selected = $bindable(),
    searchOptions = {},
    loaded = true,
  } = $props();

  // Default search configuration
  const defaultSearchOptions = $derived({
    keys: columns.map(({ key }) => key),
    threshold: 0.0,
    ignoreLocation: true,
    useExtendedSearch: true,
  });

  // Initialize Fuse.js for fuzzy search
  const fuse = $derived(
    new Fuse(data, {
      ...defaultSearchOptions,
      ...searchOptions,
    })
  );

  // State
  let query: string = $state("");

  // Filtered items based on search query
  const filteredItems = $derived(
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
