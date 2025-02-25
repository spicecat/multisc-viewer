<script lang="ts">
  import Select, { Option } from "@smui/select";

  let { groupBy = $bindable() } = $props();

</script>

<div style="margin: 1rem;">
  <div>
    <Select bind:value={groupBy} label="Group by">
      <Option value="Genotype">Genotype</Option>
      <Option value="CellType">Cell Type</Option>
    </Select>
  </div>
  <div>
    <Select bind:value={groupBy} label="Split by">
      <Option value="CellType">Genotype</Option>
      <Option value="Genotype">Cell Type</Option>
    </Select>
  </div>
</div>
