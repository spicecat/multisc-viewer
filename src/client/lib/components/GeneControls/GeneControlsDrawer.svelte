<script lang="ts">
  import GeneGroupSplit from "$lib/components/GeneControls/GeneGroupSplit.svelte";
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import Drawer, { Content } from "@smui/drawer";

  let {
    genes,
    selectedGene = $bindable(),
    loadedGenes,
    groupBy = $bindable(),
    geneControlsOpen,
  } = $props();

  const geneColumns = [{ key: "", label: "Gene" }];
</script>

<Drawer
  variant="dismissible"
  bind:open={geneControlsOpen}
>
  <Content style="display:inline-table;">
    <GeneGroupSplit bind:groupBy />
    <DataTableSearch
      label="Gene"
      data={genes}
      columns={geneColumns}
      bind:selected={selectedGene}
      loaded={loadedGenes}
    />
  </Content>
</Drawer>
<div></div>

<style lang="scss">
  // .drawer-content {
  //   display: table;
  // }

  :global(.mdc-drawer) {
    position: relative;
  }
  :global(.mdc-drawer--dismissible.mdc-drawer--open) {
    display: inline-table;
  }
</style>
