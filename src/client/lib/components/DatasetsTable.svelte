<script lang="ts">
  import DataTable, {
    Head,
    Body,
    Row,
    Cell,
    Label,
    SortValue,
  } from "@smui/data-table";
  import IconButton from "@smui/icon-button";

  let { items = [], columns = [] } = $props();
  let sort: keyof User = $state("id");
  let sortDirection: Lowercase<keyof typeof SortValue> = $state("ascending");

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
</script>

<DataTable
  sortable
  bind:sort
  bind:sortDirection
  onSMUIDataTableSorted={handleSort}
  table$aria-label="User list"
  style="width: 100%;"
  sortAscendingAriaLabel=""
  sortDescendingAriaLabel=""
>
  <Head>
    <Row>
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
    {#each items as item}
      <Row>
        {#each columns as { key }}
          <Cell numeric={typeof item[key] === "number"}>
            {item[key]}
          </Cell>
        {/each}
      </Row>
    {/each}
  </Body>
</DataTable>
