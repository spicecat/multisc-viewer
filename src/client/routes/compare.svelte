<script lang="ts">
  import ChartDisplay from "$lib/components/ChartDisplay.svelte";
  import DataTableSearch from "$lib/components/DataTableSearch.svelte";
  import meta from "$meta";
  import { dndzone } from "svelte-dnd-action";

  const { genes, gene, order }: CompareProps = $props();
  const { query } = meta;
  console.log(genes, gene, order);

  let groupBy: string = $state("Genotype"),
    selectedGene: string = $state(gene),
    dsOrder = $state(order.map((ds) => ({ id: ds })));

  const splitBy = $derived(groupBy === "Genotype" ? "CellType" : "Genotype"),
    config = $derived({ selectedGene, groupBy, splitBy });

  $effect(() => {
    if (
      !history.state ||
      history.state.selectedGene !== selectedGene ||
      history.state.groupBy !== groupBy
    )
      history.pushState(
        { selectedGene, groupBy },
        "",
        `/compare?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`
      );
  });
  const geneColumns = [{ key: "", label: "Gene" }];
</script>

<div class="row controls">
  <div class="nav">
    <a href="/" role="button">Back</a>
  </div>
  <DataTableSearch
    label="Gene"
    data={genes}
    columns={geneColumns}
    bind:selected={selectedGene}
  />
  <div class="group-split">
    <div class="row">
      <div class="col">
        <h4>Group By</h4>
        <p>{groupBy}</p>
      </div>
      <div class="col swap">
        <button onclick={() => (groupBy = splitBy)}>
          <i class="fa-solid fa-left-right"></i>
        </button>
      </div>
      <div class="col">
        <h4>Split By</h4>
        <p>{splitBy}</p>
      </div>
    </div>
  </div>
</div>
<div
  class="row charts"
  use:dndzone={{ items: dsOrder }}
  onconsider={(evt) => (dsOrder = evt.detail.items)}
  onfinalize={(evt) => (dsOrder = evt.detail.items)}
>
  {#each dsOrder as { id: dataset } (dataset)}
    <ChartDisplay {dataset} {config} />
  {/each}
</div>

<style lang="scss">
  @keyframes spin {
    0% {
      transform: rotate(0deg);
    }

    100% {
      transform: rotate(360deg);
    }
  }

  :global(.fa-spinner) {
    animation: spin 1s linear infinite;
  }

  .controls {
    height: 300px;
    gap: 1em;

    h2 {
      margin-bottom: 8px;
    }

    .gene-list {
      max-height: 180px;
      padding-block-start: 1px;
      padding-inline-start: 1px;
      padding-inline-end: 1px;
      overflow-y: scroll;

      .gene-option {
        list-style: none;
        padding: 0.2em 0.5em;
        cursor: pointer;
        opacity: 0.7;

        &:hover {
          opacity: 1;
        }

        &.selected {
          outline: 1px solid var(--primary-background);
        }
      }
    }
  }

  .charts {
    max-width: 100vw;
  }
</style>
