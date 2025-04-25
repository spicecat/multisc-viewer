<script lang="ts">
  import GeneGroupSplit from "./GeneGroupSplit.svelte";
  import DataTableSearch from "../DataTable/DataTableSearch.svelte";
  import Drawer, { Content } from "@smui/drawer";

  // Props with defaults
  let {
    genes = [],
    selectedGene = $bindable(),
    loadedGenes = true,
    groupBy = $bindable(),
    geneControlsOpen = false,
  } = $props();

  // Gene table configuration
  const geneColumns = [{ key: "", label: "Gene" }];
</script>

<Drawer
  variant="dismissible"
  bind:open={geneControlsOpen}
  class="gene-controls-drawer"
>
  <Content>
    <div class="drawer-content">
      <GeneGroupSplit bind:groupBy />
      
      <div class="gene-search">
        <DataTableSearch
          label="Gene"
          data={genes}
          columns={geneColumns}
          bind:selected={selectedGene}
          loaded={loadedGenes}
        />
      </div>
    </div>
  </Content>
</Drawer>

<style lang="scss">
  :global(.gene-controls-drawer) {
    position: relative;
    
    &.mdc-drawer--dismissible.mdc-drawer--open {
      display: inline-table;
    }
  }

  .drawer-content {
    padding: 1rem;
  }

  .gene-search {
    margin-top: 1rem;
  }
</style>
