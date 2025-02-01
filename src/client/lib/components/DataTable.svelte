<script lang="ts">
  import DataTable, { Head, Body, Row, Cell, Label } from "@smui/data-table";
  import Checkbox from "@smui/checkbox";
  import IconButton from "@smui/icon-button";

  let { items = [], columns = [] } = $props();
  let sort = $state("");
  let sortDirection: "ascending" | "descending" = $state("ascending");
  let selected = $state([]);

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
  style="width: 100%;"
  sortAscendingAriaLabel=""
  sortDescendingAriaLabel=""
>
  <Head>
    <Row>
      <Cell checkbox>
        <Checkbox />
      </Cell>
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
        <Cell checkbox>
          <Checkbox bind:group={selected} value={item} />
        </Cell>
        {#each columns as { key }}
          <Cell numeric={typeof item[key] === "number"}>
            {item[key]}
          </Cell>
        {/each}
      </Row>
    {/each}
  </Body>
</DataTable>