<script lang="ts">
  import Checkbox from "@smui/checkbox";
  import DataTable, {
    Body,
    Cell,
    Head,
    Label,
    Pagination,
    Row,
  } from "@smui/data-table";
  import IconButton from "@smui/icon-button";
  import Radio from "@smui/radio";
  import Select, { Option } from "@smui/select";

  let { items = [], columns = [], selected = $bindable() } = $props();
  let id = $derived(columns[0]?.key),
    sort = $state(""),
    sortDirection: "ascending" | "descending" = $state("ascending");

  function handleSort() {
    items.sort((a, b) => {
      const [aVal, bVal] = [a[sort], b[sort]][
        sortDirection === "ascending" ? "slice" : "reverse"
      ]();
      if (typeof aVal === "string" && typeof bVal === "string") {
        return aVal.localeCompare(bVal);
      }
      return Number(aVal) - Number(bVal);
    });
  }

  let perPage = $state(10),
    currentPage = $state(0);

  let start = $derived(currentPage * perPage),
    end = $derived(Math.min(start + perPage, items.length)),
    slice = $derived(items.slice(start, end)),
    lastPage = $derived(Math.max(Math.ceil(items.length / perPage) - 1, 0));

  $effect(() => {
    if (currentPage > lastPage) currentPage = lastPage;
  });
</script>

<DataTable
  sortable
  stickyHeader
  bind:sort
  bind:sortDirection
  onSMUIDataTableSorted={handleSort}
  table$aria-label="Data table"
  style="width: 100%;"
  sortAscendingAriaLabel=""
  sortDescendingAriaLabel=""
>
  <Head>
    <Row>
      {#if selected !== undefined}
        <Cell checkbox>
          {#if Array.isArray(selected)}
            <Checkbox />
          {/if}
        </Cell>
      {/if}
      {#each columns as { key, label, sortable = true }}
        <Cell
          columnId={key}
          numeric={typeof items[0]?.[key] === "number"}
          {sortable}
        >
          <Label>{label}</Label>
          {#if sortable && key === sort}
            <IconButton class="material-icons">
              {sortDirection === "ascending"
                ? "arrow_upward"
                : "arrow_downward"}
            </IconButton>
          {/if}
        </Cell>
      {/each}
    </Row>
  </Head>
  <Body>
    {#each slice as item (item[id])}
      <Row>
        {#if selected !== undefined}
          <Cell checkbox>
            {#if Array.isArray(selected)}
              <Checkbox bind:group={selected} value={item[id]} />
            {:else}
              <Radio bind:group={selected} value={item[id]} />
            {/if}
          </Cell>
        {/if}
        {#each columns as { key }}
          <Cell numeric={typeof item[key] === "number"}>
            {item[key]}
          </Cell>
        {/each}
      </Row>
    {/each}
  </Body>
  {#snippet paginate()}
    <Pagination>
      {#snippet rowsPerPage()}
        <Label>Rows Per Page</Label>
        <Select variant="outlined" bind:value={perPage} noLabel>
          <Option value={10}>10</Option>
          <Option value={25}>25</Option>
          <Option value={100}>100</Option>
        </Select>
      {/snippet}
      {#snippet total()}
        {start + 1}-{end} of {items.length}
      {/snippet}

      <IconButton
        class="material-icons"
        action="first-page"
        title="First page"
        onclick={() => (currentPage = 0)}
        {...{ disabled: currentPage === 0 }}
      >
        first_page
      </IconButton>
      <IconButton
        class="material-icons"
        action="prev-page"
        title="Prev page"
        onclick={() => currentPage--}
        {...{ disabled: currentPage === 0 }}
      >
        chevron_left
      </IconButton>
      <IconButton
        class="material-icons"
        action="next-page"
        title="Next page"
        onclick={() => currentPage++}
        {...{ disabled: currentPage === lastPage }}
      >
        chevron_right
      </IconButton>
      <IconButton
        class="material-icons"
        action="last-page"
        title="Last page"
        onclick={() => (currentPage = lastPage)}
        {...{ disabled: currentPage === lastPage }}
      >
        last_page
      </IconButton>
    </Pagination>
  {/snippet}
</DataTable>
