<script lang="ts">
  import TopAppBar, { Row, Section, Title } from "@smui/top-app-bar";
  import Button from "@smui/button";
  import IconButton from "@smui/icon-button";
</script>

<TopAppBar
  variant="standard"
  dense
  style="top: 0; left: 0; position: relative; margin-bottom: 1rem;"
>
  <Row>
    <Section>
      <Button href="/">
        <Title>MultiSC-Viewer</Title>
      </Button>
    </Section>
    <Section align="end">
      <Button href="/publication">
        <Title>Publications</Title>
      </Button>
      <IconButton class="material-icons" href="https://git.jasonxu.dev/JasonXu/plot-viewer/">
        source
      </IconButton>
    </Section>
  </Row>
</TopAppBar>
