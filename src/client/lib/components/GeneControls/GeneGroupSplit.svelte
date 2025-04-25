<script lang="ts">
  import Radio from "@smui/radio";
  import FormField from "@smui/form-field";

  // Props
  let { groupBy = $bindable() } = $props();

  // Available grouping options
  const groupOptions = [
    { value: "Genotype", label: "Group by Genotype" },
    { value: "CellType", label: "Group by Cell Type" }
  ];
</script>

<div class="group-controls">
  <h3 class="group-title">Grouping</h3>
  
  {#each groupOptions as { value, label }}
    <FormField>
      <Radio
        bind:group={groupBy}
        {value}
      />
      <span>{label}</span>
    </FormField>
  {/each}
</div>

<style lang="scss">
  .group-controls {
    padding: 0.5rem;
    border-radius: 4px;
    background-color: var(--surface-color, #fff);
  }

  .group-title {
    margin: 0 0 1rem;
    font-size: 1rem;
    font-weight: 500;
  }

  :global(.mdc-form-field) {
    display: block;
    margin-bottom: 0.5rem;
    
    &:last-child {
      margin-bottom: 0;
    }
  }
</style>
